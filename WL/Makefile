OPTIONS = -Wall -Wfatal-errors -O3

WL_Test.exe: WLK.o WL_Graph.o WL_Isomorph.o
	g++ -o WL_Test.exe WLK.o WL_Graph.o WL_Isomorph.o

WL_Graph.o: WL_Graph.cc WL_Graph.hh
	g++ $(OPTIONS) -c WL_Graph.cc

WL_Isomorph.o: WL_Isomorph.cc WL_Isomorph.hh WL_Graph.hh
	g++ $(OPTIONS) -c WL_Isomorph.cc

WLK.o: WLK.cpp WL_Graph.hh 
	g++ $(OPTIONS) -c WLK.cpp


clean:
	rm -f WL_Graph.o WL_Isomorph.o WL_Test.exe
