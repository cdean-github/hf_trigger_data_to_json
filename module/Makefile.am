AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64

pkginclude_HEADERS = \
  hf_trigger_data_to_json.h \
  event_display_maker.h \
  cosmic_display_maker.h

lib_LTLIBRARIES = \
  libhf_trigger_data_to_json.la

libhf_trigger_data_to_json_la_SOURCES = \
  hf_trigger_data_to_json.cc \
  event_display_maker.cc \
  cosmic_display_maker.cc

libhf_trigger_data_to_json_la_LIBADD = \
  -lKFParticle \
  -lSubsysReco \
  -lphool \
  -lcalotrigger \
  -lkfparticle_sphenix_io \
  -ltrackbase_historic

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libhf_trigger_data_to_json.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
