#! /usr/bin/env python
# encoding: utf-8

def options(opt):
	opt.recurse('hmap')
	opt.recurse('gref')
	opt.recurse('gaba')
	opt.recurse('psort')

	opt.load('compiler_c')

def configure(conf):
	conf.recurse('hmap')
	conf.recurse('gref')
	conf.recurse('gaba')
	conf.recurse('psort')

	conf.load('ar')
	conf.load('compiler_c')

	conf.env.append_value('LIB_GGSEA',
		conf.env.LIB_HMAP + conf.env.LIB_GREF + conf.env.LIB_GABA + conf.env.LIB_PSORT)
	conf.env.append_value('CFLAGS', '-O3')
	conf.env.append_value('CFLAGS', '-std=c99')
	conf.env.append_value('CFLAGS', '-march=native')


def build(bld):
	bld.recurse('hmap')
	bld.recurse('gref')
	bld.recurse('gaba')
	bld.recurse('psort')

	bld.stlib(
		source = ['ggsea.c'],
		target = 'ggsea',
		lib = bld.env.LIB_GGSEA,
		use = ['hmap', 'gref', 'gaba', 'psort'])

	bld.program(
		source = ['unittest.c'],
		target = 'unittest',
		linkflags = ['-all_load'],
		lib = bld.env.LIB_GGSEA,
		use = ['ggsea'],
		defines = ['TEST'])
