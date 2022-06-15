.PHONY: all clean

build_dir := build
target := $(build_dir)/lheureux

all: $(target)

clean:
	rm -rf $(build_dir)

$(target): src/lheureux.f
	@mkdir -p $(@D)
	gfortran -ffree-form $< -o $@

