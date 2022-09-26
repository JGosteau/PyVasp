all: split

split: procar_info.o split.o main.o
	     gcc -o split.exe procar_info.o split.o main.o	

procar_info.o: procar_info.c
	       gcc -o procar_info.o -c procar_info.c -W -Wall 

split.o: split.c procar_info.h
	 gcc -o split.o -c split.c -W -Wall 

main.o: main.c split.h
	gcc -o main.o -c main.c -W -Wall 

clean:
	rm -rf *.o

mrproper: clean
	rm -rf split.exe 
