CFLAGS=		-g -Wall -O2 -Wc++-compat #-Wextra
CPPFLAGS=	-DHAVE_KALLOC
INCLUDES=
OBJS=		kthread.o kalloc.o misc.o bseq.o sketch.o sdust.o options.o index.o \
			lchain.o align.o hit.o seed.o jump.o map.o format.o pe.o esterr.o splitidx.o \
			ksw2_ll_sse.o
PROG=		minimap2
PROG_EXTRA=	sdust minimap2-lite
LIBS=		-lm -lz -lpthread

# --- TracEon backend (optional; default off) ---
# TRACEON=1 backs minimap2's per-bucket khash_t(idx) minimizer table with
# TracEon's traceon_kmer C API (libtraceon_kmer.a, see
# https://github.com/woosflex/TracEon). All minimap2 objects still compile as
# C (gcc); only the FINAL executable link is driven by the C++ driver ($(CXX))
# so libstdc++/libtraceon_kmer symbols resolve. TRACEON_INC/TRACEON_LIB point
# at the TracEon include dir and build dir and are overridable on the command
# line. `make TRACEON=1 clean && make TRACEON=1` builds the traceon binary;
# plain `make` stays byte-for-byte on the stock khash path.
TRACEON ?= 0
CXX ?= g++
TRACEON_INC ?= $(HOME)/agent_workspace/TracEon/include
TRACEON_LIB ?= $(HOME)/agent_workspace/TracEon/build
ifeq ($(TRACEON),1)
CPPFLAGS += -DTRACEON_BACKEND -I$(TRACEON_INC)
LIBS += -L$(TRACEON_LIB) -ltraceon_kmer
endif

ifneq ($(aarch64),)
	arm_neon=1
endif

ifeq ($(arm_neon),) # if arm_neon is not defined
ifeq ($(sse2only),) # if sse2only is not defined
	OBJS+=ksw2_extz2_sse41.o ksw2_extd2_sse41.o ksw2_exts2_sse41.o ksw2_extz2_sse2.o ksw2_extd2_sse2.o ksw2_exts2_sse2.o ksw2_dispatch.o
else                # if sse2only is defined
	OBJS+=ksw2_extz2_sse.o ksw2_extd2_sse.o ksw2_exts2_sse.o
endif
else				# if arm_neon is defined
	OBJS+=ksw2_extz2_neon.o ksw2_extd2_neon.o ksw2_exts2_neon.o
    INCLUDES+=-Isse2neon
ifeq ($(aarch64),)	#if aarch64 is not defined
	CFLAGS+=-D_FILE_OFFSET_BITS=64 -mfpu=neon -fsigned-char
else				#if aarch64 is defined
	CFLAGS+=-D_FILE_OFFSET_BITS=64 -fsigned-char
endif
endif

ifneq ($(asan),)
	CFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address -ldl
endif

ifneq ($(tsan),)
	CFLAGS+=-fsanitize=thread
	LIBS+=-fsanitize=thread -ldl
endif

.PHONY:all extra clean depend
.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

extra:all $(PROG_EXTRA)

minimap2:main.o libminimap2.a
ifeq ($(TRACEON),1)
		$(CXX) $(CFLAGS) main.o -o $@ -L. -lminimap2 $(LIBS)
else
		$(CC) $(CFLAGS) main.o -o $@ -L. -lminimap2 $(LIBS)
endif

minimap2-lite:example.o libminimap2.a
ifeq ($(TRACEON),1)
		$(CXX) $(CFLAGS) $< -o $@ -L. -lminimap2 $(LIBS)
else
		$(CC) $(CFLAGS) $< -o $@ -L. -lminimap2 $(LIBS)
endif

libminimap2.a:$(OBJS)
		$(AR) -csru $@ $(OBJS)

sdust:sdust.c kalloc.o kalloc.h kdq.h kvec.h kseq.h ketopt.h sdust.h
		$(CC) -D_SDUST_MAIN $(CFLAGS) $< kalloc.o -o $@ -lz

# SSE-specific targets on x86/x86_64

ifeq ($(arm_neon),)   # if arm_neon is defined, compile this target with the default setting (i.e. no -msse2)
ksw2_ll_sse.o:ksw2_ll_sse.c ksw2.h kalloc.h
		$(CC) -c $(CFLAGS) -msse2 $(CPPFLAGS) $(INCLUDES) $< -o $@
endif

ksw2_extz2_sse41.o:ksw2_extz2_sse.c ksw2.h kalloc.h
		$(CC) -c $(CFLAGS) -msse4.1 $(CPPFLAGS) -DKSW_CPU_DISPATCH $(INCLUDES) $< -o $@

ksw2_extz2_sse2.o:ksw2_extz2_sse.c ksw2.h kalloc.h
		$(CC) -c $(CFLAGS) -msse2 -mno-sse4.1 $(CPPFLAGS) -DKSW_CPU_DISPATCH -DKSW_SSE2_ONLY $(INCLUDES) $< -o $@

ksw2_extd2_sse41.o:ksw2_extd2_sse.c ksw2.h kalloc.h
		$(CC) -c $(CFLAGS) -msse4.1 $(CPPFLAGS) -DKSW_CPU_DISPATCH $(INCLUDES) $< -o $@

ksw2_extd2_sse2.o:ksw2_extd2_sse.c ksw2.h kalloc.h
		$(CC) -c $(CFLAGS) -msse2 -mno-sse4.1 $(CPPFLAGS) -DKSW_CPU_DISPATCH -DKSW_SSE2_ONLY $(INCLUDES) $< -o $@

ksw2_exts2_sse41.o:ksw2_exts2_sse.c ksw2.h kalloc.h
		$(CC) -c $(CFLAGS) -msse4.1 $(CPPFLAGS) -DKSW_CPU_DISPATCH $(INCLUDES) $< -o $@

ksw2_exts2_sse2.o:ksw2_exts2_sse.c ksw2.h kalloc.h
		$(CC) -c $(CFLAGS) -msse2 -mno-sse4.1 $(CPPFLAGS) -DKSW_CPU_DISPATCH -DKSW_SSE2_ONLY $(INCLUDES) $< -o $@

ksw2_dispatch.o:ksw2_dispatch.c ksw2.h
		$(CC) -c $(CFLAGS) -msse4.1 $(CPPFLAGS) -DKSW_CPU_DISPATCH $(INCLUDES) $< -o $@

# NEON-specific targets on ARM

ksw2_extz2_neon.o:ksw2_extz2_sse.c ksw2.h kalloc.h
		$(CC) -c $(CFLAGS) $(CPPFLAGS) -DKSW_SSE2_ONLY -D__SSE2__ $(INCLUDES) $< -o $@

ksw2_extd2_neon.o:ksw2_extd2_sse.c ksw2.h kalloc.h
		$(CC) -c $(CFLAGS) $(CPPFLAGS) -DKSW_SSE2_ONLY -D__SSE2__ $(INCLUDES) $< -o $@

ksw2_exts2_neon.o:ksw2_exts2_sse.c ksw2.h kalloc.h
		$(CC) -c $(CFLAGS) $(CPPFLAGS) -DKSW_SSE2_ONLY -D__SSE2__ $(INCLUDES) $< -o $@

# other non-file targets

clean:
		rm -fr gmon.out *.o a.out $(PROG) $(PROG_EXTRA) *~ *.a *.dSYM build dist mappy*.so mappy.c python/mappy.c mappy.egg* .eggs

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(CPPFLAGS) -- *.c)

# DO NOT DELETE

align.o: minimap.h mmpriv.h bseq.h kseq.h ksw2.h kalloc.h
bseq.o: bseq.h kvec.h kalloc.h kseq.h
esterr.o: mmpriv.h minimap.h bseq.h kseq.h
example.o: minimap.h kseq.h
format.o: kalloc.h mmpriv.h minimap.h bseq.h kseq.h
hit.o: mmpriv.h minimap.h bseq.h kseq.h kalloc.h khash.h
index.o: kthread.h bseq.h minimap.h mmpriv.h kseq.h ksw2.h kalloc.h kvec.h
index.o: khash.h ksort.h
jump.o: mmpriv.h minimap.h bseq.h kseq.h
kalloc.o: kalloc.h
ksw2_extd2_sse.o: ksw2.h kalloc.h
ksw2_exts2_sse.o: ksw2.h kalloc.h
ksw2_extz2_sse.o: ksw2.h kalloc.h
ksw2_ll_sse.o: ksw2.h kalloc.h
kthread.o: kthread.h
lchain.o: mmpriv.h minimap.h bseq.h kseq.h kalloc.h krmq.h
main.o: bseq.h minimap.h mmpriv.h kseq.h ketopt.h
map.o: kthread.h kvec.h kalloc.h sdust.h mmpriv.h minimap.h bseq.h kseq.h
map.o: khash.h ksort.h
misc.o: mmpriv.h minimap.h bseq.h kseq.h ksort.h
options.o: mmpriv.h minimap.h bseq.h kseq.h
pe.o: mmpriv.h minimap.h bseq.h kseq.h kvec.h kalloc.h ksort.h
sdust.o: kalloc.h kdq.h kvec.h sdust.h
seed.o: mmpriv.h minimap.h bseq.h kseq.h kalloc.h ksort.h
sketch.o: kvec.h kalloc.h mmpriv.h minimap.h bseq.h kseq.h
splitidx.o: mmpriv.h minimap.h bseq.h kseq.h
