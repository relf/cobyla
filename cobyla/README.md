To generate compile commands json file required by c2rust in top directory

> mkdir build
> cd build
> cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=1 ../cobyla
> cp compile_commands.json ../cobyla

Then rust generation from C

> cd ../cobyla
> c2rust transpile compile_commands.json 