To generate compile commands json file required by c2rust in top directory

> mkdir build
> cd build
> cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=1 ../cobyla/nlopt
> cp compile_commands.json ../cobyla/nlopt/compile_commands.json

Then rust generation from C

> cd ../cobyla/nlopt
> c2rust transpile compile_commands.json 