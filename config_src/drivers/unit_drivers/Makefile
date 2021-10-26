.SILENT: # Needed to silence the directory changes

all: build/Makefile
	cd build && make
FMS: build/Makefile
	cd build && make $@
MOM%: build/Makefile
	cd build && make $@
test%: build/Makefile
	cd build && make $@

build/Makefile: CMakeLists.txt
	cmake . -Bbuild

clean:
	rm -rf build
