
all: project4 test

test: project4_test 
	./project4_test

project4_test: project4.hh rubrictest.hh project4_test.cc
	g++ -std=c++11 project4_test.cc -o project4_test

project4: project4.hh timer.hh project4_main.cc
	g++ -std=c++11 project4_main.cc -o experiment

clean:
	rm -f experiment project4_test
