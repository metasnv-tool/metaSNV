all: qaTools snpCaller
clean:
	cd src/qaTools && $(MAKE) clean
	cd src/snpCaller && $(MAKE) clean

qaTools:
	cd src/$@ && $(MAKE)
snpCaller:
	cd src/$@ && $(MAKE)

.PHONY: all clean
