@echo off
SET build_dir=build

REM 检查构建目录是否存在，如果不存在则创建它
IF NOT EXIST "%build_dir%" (
    mkdir "%build_dir%"
)

REM 切换到构建目录
cd /d "%build_dir%"

REM 运行cmake命令，这里假设CMakeLists.txt位于上级目录
cmake ..

PAUSE