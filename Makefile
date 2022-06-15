.PHONY: all clean

build_dir := build
target := $(build_dir)/lheureux
flags := -Og -fbacktrace -Wall -Wextra  \
         -fimplicit-none -g -fcheck=all -std=legacy

all: $(target)

clean:
	rm -rf $(build_dir)

$(target): src/lheureux.f
	@mkdir -p $(@D)
	gfortran -ffree-form $(flags) $< -o $@

