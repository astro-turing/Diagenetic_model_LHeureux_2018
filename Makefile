.PHONY: all clean site watch watch-pandoc watch-browser-sync

# FORTRAN PART

build_dir := build
target := $(build_dir)/lheureux
new_target := $(build_dir)/marl-pde

hdf5_cflags = $(shell h5fc -show '%' | cut -d% -f1 | cut -d' ' -f2-)
hdf5_libs = $(shell h5fc -show '%' | cut -d% -f2 | cut -d' ' -f2-)

legacy_flags := -Ofast -fbacktrace -Wall -Wextra -ffree-form \
         -fimplicit-none -g -fcheck=all -std=legacy
cflags := -Ofast -fbacktrace -Wall -Wextra  \
         -fimplicit-none -g -fcheck=all -ffree-line-length-none \
		 $(hdf5_cflags) -fintrinsic-modules-path /usr/lib64/gfortran/modules
libs := $(hdf5_libs)
object_files := $(patsubst %.f90,build/%.o,$(wildcard src/*.f90))

all: $(target) $(new_target)

build/%.o: %.f90
	@mkdir -p $(@D)
	gfortran -o $@ $(cflags) -c $<

$(target): src/lheureux.f
	@mkdir -p $(@D)
	gfortran $(legacy_flags) $< -o $@

$(new_target): $(object_files)
	@mkdir -p $(@D)
	gfortran $^ $(libs) -o $@

# ENTANGLED + PANDOC PART

theme := default
# pandoc_input := README.md $(wildcard lit/*.md)
pandoc_output := docs/index.html  docs/python-interface.html

theme_dir := .entangled/templates/$(theme)
pandoc_args += -s -t html5 -f markdown+fenced_code_attributes+fenced_divs --toc --toc-depth 2
pandoc_args += --template $(theme_dir)/template.html
pandoc_args += --css theme.css
pandoc_args += --mathjax
pandoc_args += --filter pandoc-eqnos
# pandoc_args += --syntax-definition .entangled/syntax/dhall.xml
# pandoc_args += --highlight-style $(theme_dir)/syntax.theme
pandoc_args += --section-divs
pandoc_args += --lua-filter .entangled/scripts/hide.lua
pandoc_args += --lua-filter .entangled/scripts/annotate.lua
pandoc_args += --lua-filter .entangled/scripts/make.lua
static_files := $(theme_dir)/theme.css $(theme_dir)/static
static_targets := $(static_files:$(theme_dir)/%=docs/%)
functional_deps := Makefile $(wildcard .entangled/scripts/*.lua) $(theme_dir)/template.html $(theme_dir)/syntax.theme

# figure_src := $(wildcard fig/*)
# figure_targets := $(figure_src:%=docs/%)

site: $(pandoc_output) $(static_targets) $(figure_targets)

clean:
	rm -rf docs $(build_dir)

# $(figure_targets): docs/fig/%: fig/%
# 	@mkdir -p $(@D)
# 	cp $< $@

$(static_targets): docs/%: $(theme_dir)/%
	@mkdir -p $(@D)
	rm -rf $@
	cp -r $< $@

docs/%.html: lit/%.md $(functional_deps)
	@mkdir -p $(@D)
	pandoc $(pandoc_args) -o $@ $<

# Starts a tmux with Entangled, Browser-sync and an Inotify loop for running
# Pandoc.
watch:
	@tmux new-session make --no-print-directory watch-pandoc \; \
		split-window -v make --no-print-directory watch-browser-sync \; \
		split-window -v entangled daemon \; \
		select-layout even-vertical \;

watch-pandoc:
	@while true; do \
		inotifywait -e close_write -r .entangled Makefile README.md lit; \
		make site; \
	done

watch-browser-sync:
	browser-sync start -w -s docs
