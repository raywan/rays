@echo off
if not exist build mkdir build
pushd build
REM NOTE(ray): This will try to compile the non-Windows main file
cl /O2 /EHsc /Zi /I ..\include /c ..\src\*.cpp
REM NOTE(ray): Just manually type in each object file for now
link win32_main.obj render.obj mesh.obj world.obj utils.obj metrics.obj ray.obj bvh.obj camera.obj primitive.obj /OUT:rays.exe
popd
