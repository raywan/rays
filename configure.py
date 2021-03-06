#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import sys
import fnmatch
from ninja_syntax import Writer

PROJECT_NAME = 'toy' if sys.platform != 'win32' else 'toy.exe'

CC = 'g++' if sys.platform != 'win32' else 'cl'
CFLAGS = []
if sys.platform == 'win32':
    CFLAGS = ['-O2', '-EHsc', '-Zo', '/fp:fast', '-Iinclude']
else:
    CFLAGS = ['-std=c++11', '-O3', '-pthread', '-march=native', '-Iinclude']

SRC_DIR = 'src'
BUILD_DIR = 'build'
BIN_DIR = 'bin'

with open('build.ninja', 'w') as build_file:
    n = Writer(build_file)

    n.comment('THIS FILE IS GENERATED BY configure.py')
    n.comment('EDITS WILL BE OVERWRITTEN')
    n.newline()

    ############################################################################
    # VARIABLES
    ############################################################################

    n.variable(key='ninja_required_version', value='1.9')

    if sys.platform == 'win32':
        n.variable(key='msvc_deps_prefix', value='Note: including file:')

    n.variable(key='cc', value=CC)
    n.variable(key='cflags', value=' '.join(CFLAGS))
    n.variable(key='project_name', value=PROJECT_NAME)
    n.variable(key='src_dir', value=SRC_DIR)
    n.variable(key='bin_dir', value=BIN_DIR)
    n.variable(key='build_dir', value=BUILD_DIR)

    n.newline()

    ############################################################################
    # RULES
    ############################################################################

    if sys.platform == 'win32':
        n.rule('compile', command='$cc /showIncludes $cflags -c $in -Fo$out', deps='msvc')
    else:
        n.rule('compile',
               command='$cc $cflags -c $in -o $out -MMD -MF $out.d',
               depfile='$out.d')

    if sys.platform == 'win32':
        n.rule('link', command='link $in /OUT:$out')
    else:
        n.rule('link', command='$cc $in -o $bin_dir/$out')

    n.rule('clean_all', command='rm -rf $build_dir $bin_dir $project_name')
    if sys.platform == 'win32':
        n.rule('make_dirs', command='cmd /c if not exist $build_dir mkdir $build_dir')
    else:
        n.rule('make_dirs', command='mkdir -p $build_dir $bin_dir')
    if sys.platform == 'win32':
        n.rule('create_sym_link', command='cmd /c mklink $project_name.exe $build_dir\$project_name.exe')
    else:
        n.rule('create_sym_link', command='ln -sf $bin_dir/$project_name $project_name')

    n.newline()

    ############################################################################
    # BUILDS
    ############################################################################

    n.build(outputs='dirs', rule='make_dirs')
    n.build(outputs='fresh', rule='clean_all')
    n.build(outputs='sym', rule='create_sym_link')

    sources = []
    for (root, dirnames, filenames) in os.walk(SRC_DIR):
        for filename in fnmatch.filter(filenames, '*.cpp'):
            if sys.platform != 'win32' and 'win32_' not in filename:
                sources.append(os.path.join(root, filename))
            elif sys.platform == 'win32' and filename != 'main.cpp':
                sources.append(os.path.join(root, filename))


    for source in sources:
        if sys.platform == 'win32':
            n.build(outputs=source.replace('.cpp', '.obj').replace('src',
                    BUILD_DIR), rule='compile', inputs=source)
        else:
            n.build(outputs=source.replace('.cpp', '.o').replace('src',
                    BUILD_DIR), rule='compile', inputs=source)


    o_files = []
    if sys.platform == 'win32':
        o_files = [x.replace('.cpp', '.obj').replace('src', BUILD_DIR) for x in sources]
    else:
        o_files = [x.replace('.cpp', '.o').replace('src', BUILD_DIR) for x in sources]

    n.build(outputs=PROJECT_NAME, rule='link', inputs=o_files)

    n.newline()

    ############################################################################
    # DEFAULTS
    ############################################################################

    n.default('dirs')
    n.default(PROJECT_NAME)
    if sys.platform != 'win32':
        n.default('sym')
