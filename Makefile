simple_example: obj/simple_example.o obj/fuzzy.o
	gcc -o bin/simple_example obj/simple_example.o obj/fuzzy.o -lm

simple_example2: obj/simple_example2.o obj/fuzzy.o
	gcc -o bin/simple_example2 obj/simple_example2.o obj/fuzzy.o -lm

obj/simple_example.o: src/simple_example.c include/fuzzy.h
	gcc -c src/simple_example.c -I include

obj/fuzzy.o: src/fuzzy.c include/fuzzy.h
	gcc -c fuzzy.c -I include

clean:
	rm bin/simple_example bin/simple_example2\
		obj/simple_example.o obj/simple_example2\
		obj/fuzzy.o
