A *very* stripped down library for fuzzy logic. This is planned to be the core of
a genetic fuzzy library, so only basic needs are supplied. Only supports triangular
mambership functions and "mean of moments" defuzzification. Emphasis is put on speed
and simple interface (though care was taken to also be correct and safe in the process).

Makefile is planned, but not here yet. Compile with
```gcc -o fuzzy.o -c fuzzy.c -lm```
```gcc -o prog program.c fuzzy.o -lm```
