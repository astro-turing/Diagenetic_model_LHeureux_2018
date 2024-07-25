build_dir := build
target := $(build_dir)/marlpde08

f08_flags = -O3 -ffree-form -std=f2008 -g

$(target): src/lheureux.f90
	@mkdir -p $(@D)
	gfortran $(f08_flags) $< -o $@

.PHONY: clean 

clean :
	rm -rf $(build_dir) *.mod