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
    if env['PLATFORM']=='win32' and env.get('CONTAINERNAME','')=='GlastRelease':
        env.Tool('findPkgPath', package = 'FluxSvc') 
def exists(env):
    return 1;
