hmm_dir = /scratch/fishjord/apps/hmmer_svn/
#hmm_dir = /scratch/fishjord/apps/hmmer_svn_linux/
jvm_home = /System/Library/Java/JavaVirtualMachines/1.6.0.jdk/Contents/Home/
#jvm_home = /usr/java/latest
inc_dirs = -I$(hmm_dir)/easel/ -I$(hmm_dir)/src/ -I$(jvm_home)/include -I$(jvm_home)/include/linux
lib_dirs = -L$(hmm_dir)/easel/ -L$(hmm_dir)/src/
libs = -lhmmer -leasel -lm

gcc_opts = -std=gnu99 -O3 -fomit-frame-pointer -fstrict-aliasing -mpentiumpro -msse2 #-ggdb

libhmmerwrapper.dylib: hmmer_wrapper.c edu_msu_cme_rdp_alignment_hmm_jni_HMMER3.c
	gcc -dynamiclib -o libhmmerwrapper.dylib \
	$(gcc_opts) \
	edu_msu_cme_rdp_alignment_hmm_jni_HMMER3.c hmmer_wrapper.c \
	$(inc_dirs) $(lib_dirs) $(libs)

libhmmerwrapper.so: libhmmer.so libeasel.so hmmer_wrapper.c edu_msu_cme_rdp_alignment_hmm_jni_HMMER3.c
	gcc -shared -Wl,-soname,libhmmerwrapper.so -o libhmmerwrapper.so \
	$(gcc_opts) -fPIC \
	edu_msu_cme_rdp_alignment_hmm_jni_HMMER3.c hmmer_wrapper.c \
	$(inc_dirs) -L. $(libs) -lc

libhmmer.so: 
	gcc -shared -Wl,-soname,libhmmer.so -o libhmmer.so $(hmm_dir)/src/*.o $(hmm_dir)/src/impl/*.o
libeasel.so: 
	gcc -Wl,-soname,libeasel.so -o libeasel.so -shared $(hmm_dir)/easel/*.o

wrapper: hmmer_wrapper.c hmmer_wrapper.h
	gcc -o wrapper $(gcc_opts) \
	hmmer_wrapper.c \
	$(inc_dirs) $(lib_dirs) $(libs)
