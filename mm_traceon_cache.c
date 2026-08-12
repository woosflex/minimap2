/*
 * mm_traceon_cache.c — TRC1 ".tcache" flat-array index cache (TracEon backend)
 *
 * OPT-IN, additive cache mode for RECURRING runs: alongside the stock khash
 * backend (plain `make`) and the traceon table backend (`make TRACEON=1`),
 * this file implements a third, mmap-loadable flat format that eliminates the
 * load-time table rebuild entirely. Both khash and traceon REBUILD their
 * per-bucket tables on load (every .mmi record is inserted into a hash table);
 * a .tcache file instead stores each bucket's minimizer entries as a SORTED
 * flat (key,value) array plus cumulative offset tables, so loading is: mmap +
 * one whole-file CRC32C check + pointer fixups. No inserts, no table rebuild,
 * no khash/traceon objects. Lookups use binary search inside the bucket
 * (O(log n) per minimizer) and reproduce the exact khash semantics
 * (idx_hash/idx_eq ignore key bit 0 — the "singleton" flag — so a lookup
 * with bit0=0 finds an entry stored with bit0=1; values are packed exactly as
 * khash stores them: singleton -> position u64, multi -> (start_p<<32)|n).
 *
 * The format is only compiled into TRACEON builds (same TRACEON_BACKEND
 * guard as the traceon table backend); the stock build never sees it.
 *
 * ---------------------------------------------------------------------------
 * TRC1 FILE FORMAT (little-endian, host byte order — same convention as .mmi)
 * ---------------------------------------------------------------------------
 *   Header (72 bytes):
 *     u32  magic     "TRC1"
 *     u16  version   1
 *     u16  flags     bit0 has_seq, bit1 has_names (mirrors the stored MM_I_* flags)
 *     i32  w         minimizer window
 *     i32  k         minimizer k-mer size
 *     i32  b         bucket bits (n_bucket == 1<<b)
 *     u32  n_seq     number of reference sequences
 *     u32  n_bucket  1<<b
 *     u64  sum_len   total reference length in bases
 *     u64  n_entry   total (key,value) entries across ALL buckets
 *     u64  names_len total bytes of the names blob (incl. per-name length bytes)
 *     u64  p_total   total u64 count across ALL per-bucket position (p) arrays
 *     u32  idx_flag  the MM_I_* flags the index was built with
 *     u32  reserved  0
 *     u32  pad       0
 *
 *   Reference block (restores mi->seq and mi->S WITHOUT parsing FASTA):
 *     names blob  names_len bytes: per sequence [u8 len][len bytes] (no NUL)
 *     (pad to 8)
 *     lens array  n_seq * u32 sequence lengths (offsets = cumulative sum)
 *     S blob      has_seq ? ((sum_len+7)/8)*4 bytes : 0 — the 4-bit packed
 *                 reference, byte-identical layout to mm_idx_t::S
 *     (pad to 8)
 *
 *   Table block (zero-rebuild, mmap-pointable):
 *     p_off   (n_bucket+1) * u64 cumulative counts of p-array entries
 *             (p_off[0]=0; bucket i's p array = p_blob[p_off[i]..p_off[i+1]))
 *     e_off   (n_bucket+1) * u64 cumulative counts of (key,value) entries
 *             (bucket i's entries = entries[e_off[i]..e_off[i+1]))
 *     entries n_entry * (u64 key, u64 value), per-bucket CONTIGUOUS and SORTED
 *             by key>>1 (raw key order == key>>1 order, since every bucket has
 *             at most one entry per key>>1). Binary-searchable per bucket.
 *     p_blob  p_total * u64 concatenated per-bucket position arrays (each
 *             bucket's is already sorted by position, as in the khash build)
 *
 *   Trailer:
 *     u32  crc32c  CRC-32C (Castagnoli, init 0xFFFFFFFF, final XOR) over the
 *             WHOLE file from byte 0 up to (excluding) this trailer — the
 *             TracEon .traceon v4 whole-payload integrity pattern. Computed
 *             with TracEon's include/Crc32c.h (see mm_traceon_crc.cpp).
 *
 * The layout is computed by tcache_calc_layout() — the ONLY place the offset
 * math lives — used by both save and load so they can never drift apart.
 *
 * ---------------------------------------------------------------------------
 * USAGE
 * ---------------------------------------------------------------------------
 *   ./minimap2 -d ref.tcache ref.fa          # one-time build (slow path)
 *   ./minimap2 ref.tcache reads.fq           # recurring runs: ~zero-rebuild
 *
 * `.tcache` extension on -d selects the format; a TRC1 magic on input selects
 * the loader. Both are dispatched from index.c (mm_idx_reader_*).
 */

#ifdef TRACEON_BACKEND

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <assert.h>

#include "minimap.h"
#include "mmpriv.h"
#include "ksort.h"
#include "kalloc.h" // kmalloc/kcalloc (NULL km falls back to malloc/calloc)
#include "kmerindex_c_api.h"

/* ---- CRC-32C (Castagnoli) — C shim over TracEon's include/Crc32c.h, which is
 * ---- header-only C++; the thin wrapper lives in mm_traceon_crc.cpp and is
 * ---- linked via $(CXX) in TRACEON builds (same link model as libtraceon_kmer). */
typedef struct mm_crc32c_s mm_crc32c_s;
mm_crc32c_s *mm_crc32c_new(void);
void mm_crc32c_free(mm_crc32c_s *s);
void mm_crc32c_update(mm_crc32c_s *s, const void *data, size_t len);
uint32_t mm_crc32c_final(mm_crc32c_s *s);
uint32_t mm_crc32c(const void *data, size_t len);

#define TCACHE_MAGIC     "TRC1"
#define TCACHE_VERSION   1
#define TCACHE_HDR_SIZE  72

#define TCACHE_F_HAS_SEQ   0x1
#define TCACHE_F_HAS_NAMES 0x2

#define TC_OFF_MAGIC    0
#define TC_OFF_VERSION  4
#define TC_OFF_FLAGS    6
#define TC_OFF_W        8
#define TC_OFF_K        12
#define TC_OFF_B        16
#define TC_OFF_NSEQ     20
#define TC_OFF_NBUCKET  24
#define TC_OFF_SUMLEN   28
#define TC_OFF_NENTRY   36
#define TC_OFF_NAMESLEN 44
#define TC_OFF_PTOTAL   52
#define TC_OFF_IDXFLAG  60
#define TC_OFF_RSV      64
#define TC_OFF_PAD      68

/* Byte layout of everything after the header (offsets relative to file base). */
typedef struct {
	uint64_t names_off;  // names blob (u8-len-prefixed)
	uint64_t lens_off;   // n_seq * u32 lengths (8-aligned)
	uint64_t s_off;      // packed reference (S) blob
	uint64_t s_size;     // bytes of S blob
	uint64_t poff_off;   // (n_bucket+1) * u64 cumulative p counts
	uint64_t eoff_off;   // (n_bucket+1) * u64 cumulative entry counts
	uint64_t ent_off;    // n_entry * 16 bytes (key,value)
	uint64_t pblob_off;  // p_total * 8 bytes
	uint64_t payload_size; // == file size minus the 4-byte CRC trailer
} tcache_layout_t;

static void tcache_calc_layout(uint64_t sum_len, uint32_t has_seq, uint32_t n_seq,
                               uint64_t names_len, uint32_t n_bucket,
                               uint64_t n_entry, uint64_t p_total,
                               tcache_layout_t *L)
{
	L->names_off = TCACHE_HDR_SIZE;
	L->lens_off = (L->names_off + names_len + 7) & ~(uint64_t)7;
	L->s_off = L->lens_off + (uint64_t)n_seq * 4;
	L->s_size = has_seq? ((sum_len + 7) / 8) * 4 : 0;
	L->poff_off = (L->s_off + L->s_size + 7) & ~(uint64_t)7;
	L->eoff_off = L->poff_off + (uint64_t)(n_bucket + 1) * 8;
	L->ent_off = L->eoff_off + (uint64_t)(n_bucket + 1) * 8;
	L->pblob_off = L->ent_off + n_entry * 16;
	L->payload_size = L->pblob_off + p_total * 8;
}

static uint16_t rd16(const uint8_t *p) { uint16_t v; memcpy(&v, p, 2); return v; }
static uint32_t rd32(const uint8_t *p) { uint32_t v; memcpy(&v, p, 4); return v; }
static uint64_t rd64(const uint8_t *p) { uint64_t v; memcpy(&v, p, 8); return v; }
static void wr16(uint8_t *p, uint16_t v) { memcpy(p, &v, 2); }
static void wr32(uint8_t *p, uint32_t v) { memcpy(p, &v, 4); }
static void wr64(uint8_t *p, uint64_t v) { memcpy(p, &v, 8); }

/* Per-bucket entry collection + radix sort by raw key (== key>>1 order). */
typedef struct { uint64_t k, v; } tcache_pair_t;
#define sort_key_tcp(a) ((a).k)
KRADIX_SORT_INIT(tcp, tcache_pair_t, sort_key_tcp, 8)

/* Write `n` bytes from `buf` and feed them to the running CRC. */
static int tcache_write(FILE *fp, mm_crc32c_s *crc, const void *buf, size_t n)
{
	if (fwrite(buf, 1, n, fp) != n) return -1;
	mm_crc32c_update(crc, buf, n);
	return 0;
}

int mm_tcache_dump(FILE *fp, const mm_idx_t *mi)
{
	uint64_t i, j, sum_len = 0, names_len = 0, n_entry = 0, p_total = 0;
	uint32_t has_seq = !(mi->flag & MM_I_NO_SEQ);
	uint32_t has_names = !(mi->flag & MM_I_NO_NAME);
	uint32_t n_bucket = 1U << mi->b;
	uint64_t *ecount, *pcount;
	tcache_layout_t L;
	mm_crc32c_s *crc = 0;
	uint8_t hdr[TCACHE_HDR_SIZE];

	if (fp == 0 || mi == 0) return -1;

	/* pass A: sizes */
	for (i = 0; i < mi->n_seq; ++i) {
		sum_len += mi->seq[i].len;
		if (has_names && mi->seq[i].name)
			names_len += 1 + strlen(mi->seq[i].name);
	}
	ecount = (uint64_t*)calloc(n_bucket, 8);
	pcount = (uint64_t*)calloc(n_bucket, 8);
	if (ecount == 0 || pcount == 0) goto fail0;
	for (i = 0; i < n_bucket; ++i) {
		const mm_idx_bucket_t *b = &mi->B[i];
		ecount[i] = mi->is_tcache? (uint64_t)b->ne
		                        : (b->h? kmerindex_size((const kmerindex_t*)b->h) : 0);
		pcount[i] = (uint64_t)b->n;
		n_entry += ecount[i];
		p_total += pcount[i];
	}
	tcache_calc_layout(sum_len, has_seq, mi->n_seq, names_len, n_bucket, n_entry, p_total, &L);

	/* header */
	memset(hdr, 0, sizeof(hdr));
	memcpy(hdr + TC_OFF_MAGIC, TCACHE_MAGIC, 4);
	wr16(hdr + TC_OFF_VERSION, TCACHE_VERSION);
	wr16(hdr + TC_OFF_FLAGS, (has_seq? TCACHE_F_HAS_SEQ : 0) | (has_names? TCACHE_F_HAS_NAMES : 0));
	wr32(hdr + TC_OFF_W, (uint32_t)mi->w);
	wr32(hdr + TC_OFF_K, (uint32_t)mi->k);
	wr32(hdr + TC_OFF_B, (uint32_t)mi->b);
	wr32(hdr + TC_OFF_NSEQ, mi->n_seq);
	wr32(hdr + TC_OFF_NBUCKET, n_bucket);
	wr64(hdr + TC_OFF_SUMLEN, sum_len);
	wr64(hdr + TC_OFF_NENTRY, n_entry);
	wr64(hdr + TC_OFF_NAMESLEN, names_len);
	wr64(hdr + TC_OFF_PTOTAL, p_total);
	wr32(hdr + TC_OFF_IDXFLAG, (uint32_t)mi->flag);

	crc = mm_crc32c_new();
	if (crc == 0) goto fail0;
	if (tcache_write(fp, crc, hdr, sizeof(hdr)) != 0) goto fail;

	/* names blob: [u8 len][len bytes] per sequence, no NULs (same as .mmi) */
	for (i = 0; i < mi->n_seq; ++i) {
		if (has_names && mi->seq[i].name) {
			uint8_t nb[256];
			size_t l = strlen(mi->seq[i].name);
			assert(l < 256); // same u8-len constraint as the .mmi format
			nb[0] = (uint8_t)l;
			memcpy(nb + 1, mi->seq[i].name, l);
			if (tcache_write(fp, crc, nb, l + 1) != 0) goto fail;
		} else {
			uint8_t z = 0;
			if (tcache_write(fp, crc, &z, 1) != 0) goto fail;
		}
	}
	/* pad names blob to 8 */
	{
		uint8_t zb[8] = {0,0,0,0,0,0,0,0};
		uint64_t pad = L.lens_off - (TCACHE_HDR_SIZE + names_len);
		if (pad && tcache_write(fp, crc, zb, (size_t)pad) != 0) goto fail;
	}
	/* lens array */
	if (mi->n_seq) {
		uint32_t *lens = (uint32_t*)malloc(mi->n_seq * 4);
		if (lens == 0) goto fail;
		for (i = 0; i < mi->n_seq; ++i) lens[i] = mi->seq[i].len;
		if (tcache_write(fp, crc, lens, mi->n_seq * 4) != 0) { free(lens); goto fail; }
		free(lens);
	}
	/* packed reference (S) — identical byte layout to mm_idx_t::S */
	if (has_seq) {
		size_t n32 = (size_t)((sum_len + 7) / 8);
		if (fwrite(mi->S, 4, n32, fp) != n32) goto fail;
		mm_crc32c_update(crc, mi->S, n32 * 4);
	}
	/* pad S blob to 8 (before the u64 tables) */
	{
		uint8_t zb[8] = {0,0,0,0,0,0,0,0};
		uint64_t pad = L.poff_off - (L.s_off + L.s_size);
		if (pad && tcache_write(fp, crc, zb, (size_t)pad) != 0) goto fail;
	}
	/* p_off table: cumulative p-array entry counts */
	for (i = 0, j = 0; i <= n_bucket; ++i) {
		uint64_t v = j;
		if (tcache_write(fp, crc, &v, 8) != 0) goto fail;
		if (i < n_bucket) j += pcount[i];
	}
	/* e_off table: cumulative (key,value) entry counts */
	for (i = 0, j = 0; i <= n_bucket; ++i) {
		uint64_t v = j;
		if (tcache_write(fp, crc, &v, 8) != 0) goto fail;
		if (i < n_bucket) j += ecount[i];
	}
	/* entries blob: per bucket, contiguous and sorted by key */
	for (i = 0; i < n_bucket; ++i) {
		uint64_t c;
		tcache_pair_t *pairs;
		if (ecount[i] == 0) continue;
		pairs = (tcache_pair_t*)malloc(ecount[i] * sizeof(tcache_pair_t));
		if (pairs == 0) goto fail;
		if (mi->is_tcache) {
			// flat source: already sorted by construction
			for (c = 0; c < ecount[i]; ++c) {
				pairs[c].k = mi->B[i].fe[c<<1];
				pairs[c].v = mi->B[i].fe[(c<<1)+1];
			}
		} else {
			kmerindex_iter_t it;
			uint64_t kk, vv;
			kmerindex_iter_begin((const kmerindex_t*)mi->B[i].h, &it);
			for (c = 0; kmerindex_iter_next(&it, &kk, &vv); ++c) {
				pairs[c].k = kk;
				pairs[c].v = vv;
			}
			assert(c == ecount[i]);
			radix_sort_tcp(pairs, pairs + c);
		}
		if (fwrite(pairs, 16, ecount[i], fp) != ecount[i]) { free(pairs); goto fail; }
		mm_crc32c_update(crc, pairs, ecount[i] * 16);
		free(pairs);
	}
	/* p_blob: concatenated per-bucket position arrays (already sorted) */
	for (i = 0; i < n_bucket; ++i) {
		if (pcount[i] == 0) continue;
		if (fwrite(mi->B[i].p, 8, pcount[i], fp) != pcount[i]) goto fail;
		mm_crc32c_update(crc, mi->B[i].p, pcount[i] * 8);
	}
	/* whole-file CRC32C trailer */
	{
		uint32_t c = mm_crc32c_final(crc);
		if (fwrite(&c, 4, 1, fp) != 1) goto fail;
	}
	fflush(fp);
	mm_crc32c_free(crc);
	free(ecount); free(pcount);
	return 0;

fail:
	mm_crc32c_free(crc);
fail0:
	free(ecount); free(pcount);
	return -1;
}

mm_idx_t *mm_tcache_load(FILE *fp)
{
	int fd;
	struct stat st;
	uint8_t *map;
	uint64_t size;
	tcache_layout_t L;
	mm_idx_t *mi;
	uint32_t n_bucket, has_seq, has_names, idx_flag;
	uint64_t sum_len, n_entry, names_len, p_total;
	uint32_t w, k, b, n_seq, i;
	double t0, t_mmap, t_crc;

	if (fp == 0) return 0;
	fd = fileno(fp);
	if (fd < 0 || fstat(fd, &st) != 0) return 0;
	size = (uint64_t)st.st_size;
	if (size < TCACHE_HDR_SIZE + 4) return 0;
#ifdef WIN32
	if (_ftelli64(fp) >= (int64_t)size) return 0; // stream already consumed
#else
	if (ftello(fp) >= (off_t)size) return 0;      // stream already consumed (multi-part EOF, like mm_idx_load)
#endif
	t0 = realtime();
	map = (uint8_t*)mmap(0, size, PROT_READ, MAP_PRIVATE, fd, 0);
	if (map == MAP_FAILED) return 0;
	t_mmap = realtime();

	/* magic + version */
	if (memcmp(map + TC_OFF_MAGIC, TCACHE_MAGIC, 4) != 0) goto fail;
	if (rd16(map + TC_OFF_VERSION) != TCACHE_VERSION) goto fail;

	/* whole-file CRC32C trailer: covers [0, size-4) — TracEon v4 integrity pattern */
	{
		uint32_t want = rd32(map + size - 4);
		uint32_t got = mm_crc32c(map, size - 4);
		if (want != got) {
			fprintf(stderr, "[ERROR] mm_tcache_load: CRC32C mismatch (%s): stored %08x, computed %08x\n",
				"file corrupted or truncated", want, got);
			goto fail;
		}
	}
	t_crc = realtime();

	w = rd32(map + TC_OFF_W); k = rd32(map + TC_OFF_K); b = rd32(map + TC_OFF_B);
	n_seq = rd32(map + TC_OFF_NSEQ); n_bucket = rd32(map + TC_OFF_NBUCKET);
	sum_len = rd64(map + TC_OFF_SUMLEN); n_entry = rd64(map + TC_OFF_NENTRY);
	names_len = rd64(map + TC_OFF_NAMESLEN); p_total = rd64(map + TC_OFF_PTOTAL);
	idx_flag = rd32(map + TC_OFF_IDXFLAG);
	has_seq = (rd16(map + TC_OFF_FLAGS) & TCACHE_F_HAS_SEQ) != 0;
	has_names = (rd16(map + TC_OFF_FLAGS) & TCACHE_F_HAS_NAMES) != 0;

	if (b < 1 || b > 30) goto fail;                 // avoid 1U<<b overflow
	if (n_bucket != 1U << b) goto fail;
	if (n_seq != 0 && sum_len == 0) goto fail;      // sanity
	tcache_calc_layout(sum_len, has_seq, n_seq, names_len, n_bucket, n_entry, p_total, &L);
	if (L.payload_size + 4 != size) goto fail;      // exact layout; CRC already guards content

	mi = mm_idx_init(w, k, b, idx_flag);
	if (mi == 0) goto fail;
	mi->is_tcache = 1;
	mi->tcache_map = map;
	mi->tcache_size = (int64_t)size;
	mi->n_seq = n_seq;
	mi->seq = (mm_idx_seq_t*)kcalloc(mi->km, n_seq, sizeof(mm_idx_seq_t));

	/* reference block: names + lengths + packed sequence (NO FASTA parsing) */
	{
		const uint8_t *np = map + L.names_off;
		const uint8_t *lp = map + L.lens_off;
		uint64_t sum = 0;
		for (i = 0; i < n_seq; ++i) {
			mm_idx_seq_t *s = &mi->seq[i];
			uint8_t l = np[0];
			if (has_names && l) {
				s->name = (char*)kmalloc(mi->km, l + 1);
				memcpy(s->name, np + 1, l);
				s->name[l] = 0;
			} else s->name = 0;
			np += 1 + l;
			s->len = rd32(lp + (uint64_t)i * 4);
			s->offset = sum;
			s->is_alt = 0;
			sum += s->len;
		}
	}
	if (has_seq) mi->S = (uint32_t*)(map + L.s_off);

	/* table block: point every bucket at its mmapped arrays — ZERO inserts */
	{
		const uint64_t *poff = (const uint64_t*)(map + L.poff_off);
		const uint64_t *eoff = (const uint64_t*)(map + L.eoff_off);
		const uint64_t *ent  = (const uint64_t*)(map + L.ent_off);
		const uint64_t *pblob= (const uint64_t*)(map + L.pblob_off);
		for (i = 0; i < n_bucket; ++i) {
			mm_idx_bucket_t *b = &mi->B[i];
			uint64_t pe = poff[i+1] - poff[i];
			uint64_t ee = eoff[i+1] - eoff[i];
			b->n = (int32_t)pe;
			if (pe) b->p = (uint64_t*)(pblob + poff[i]);
			if (ee) { b->fe = ent + (eoff[i]<<1); b->ne = (int32_t)ee; }
		}
	}

	/* advance the FILE position to EOF so mm_idx_reader_eof() terminates the
	 * multi-part loop exactly like a consumed .mmi stream */
#ifdef WIN32
	_fseeki64(fp, (int64_t)size, SEEK_SET);
#else
	fseeko(fp, (off_t)size, SEEK_SET);
#endif
	if (mm_verbose >= 3) {
		double t1 = realtime();
		fprintf(stderr, "[M::%s] mmap %.1f ms, crc32c %.1f ms, fixup %.1f ms (total %.1f ms) — %lld entries, %lld p-entries, %lld ref bases\n",
			__func__, (t_mmap - t0) * 1e3, (t_crc - t_mmap) * 1e3, (t1 - t_crc) * 1e3,
			(t1 - t0) * 1e3, (long long)n_entry, (long long)p_total, (long long)sum_len);
	}
	return mi;

fail:
	munmap(map, size);
	return 0;
}

#endif /* TRACEON_BACKEND */
