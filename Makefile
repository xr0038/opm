all: src example

src: always
	make -C src
	@if [ ! -d include ]; then mkdir ./include; fi
	rsync -avu src/*.h ./include
	@if [ ! -d lib ]; then mkdir ./lib; fi
	rsync -avu src/*.a ./lib

example: always
	make -C example


always:
