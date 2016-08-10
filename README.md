A *very* stripped down library for fuzzy logic. This is planned to be the core of
a genetic fuzzy library, so only basic needs are supplied. Only supports triangular
mambership functions and "mean of moments" defuzzification. Emphasis is put on speed
and simple interface (though care was taken to also be correct and safe in the process).

GNU Autotools were used to create project Makefiles. Tested on Ubuntu 14.04 and FreeBSD 10.4.
For some reason, in *BSD, `configure` is not made executable, so `chmod +x configure` is necessary.
Autoconf supports parallel builds, and it is recomended to `mkdir build && cd build && ../configure`
to build and then test the library.


