# Based on Johan's makefile 

build_dir := build
target := $(build_dir)/lheureux


legacy_flags := -Ofast -fbacktrace -Wall -Wextra -ffree-form \
         -fimplicit-none -g -fcheck=all -std=legacy -Wno-tabs
cflags := -Ofast -fbacktrace -Wall -Wextra  \
         -fimplicit-none -g -fcheck=all -ffree-line-length-none \
		  -fintrinsic-modules-path /usr/lib64/gfortran/modules

object_files := $(patsubst %.f,build/%.o,$(wildcard src/*.f))

all: $(target) $(new_target)

build/%.o: %.f
	@mkdir -p $(@D)
	gfortran -o $@ $(cflags) -c $<

$(target): src/lheureux.f
	@mkdir -p $(@D)
	gfortran $(legacy_flags) $< -o $@

$(new_target): $(object_files)
	@mkdir -p $(@D)
	gfortran $(libs) $^ -o $@

