.PHONY: test docs

test:
	julia test/runtests.jl

docs:
	julia --project=docs/ docs/make.jl
