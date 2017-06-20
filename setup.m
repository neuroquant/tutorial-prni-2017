HOME=getenv('HOME');
USE_OCTAVE=1;

basedir=fullfile(HOME,'prni2017-site-effects');
addpath(fullfile(basedir,'data')); 
addpath(fullfile(basedir,'src')); 
addpath(genpath(fullfile(basedir,'src','external')));
addpath(fullfile(basedir,'src','mwwTest'));

if(USE_OCTAVE)
    setup_octave;
end