# $Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['CRflux'])
    env.Tool('addLibrary', library = env['cfitsioLibs'])
    env.Tool('addLibrary', library = env['clhepLibs'])
    env.Tool('EventLib')
    env.Tool('addLibrary', library = env['xercesLibs'])
    env.Tool('astroLib')
    env.Tool('xmlBaseLib')
    env.Tool('fluxLib')
def exists(env):
    return 1;
