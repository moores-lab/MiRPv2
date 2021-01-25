#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import stat


#make mirp files executable
os.chmod('mirp/generate_seam_references.py', stat.S_IRWXU)
os.chmod('mirp/generate_segment_averages.py', stat.S_IRWXU)
os.chmod('mirp/mirp_pf_sorting', stat.S_IRWXU)
os.chmod('mirp/mirp_initial_seam', stat.S_IRWXU)
os.chmod('mirp/mirp_seam_check', stat.S_IRWXU)
os.chmod('mirp/plot_eulerxy.py', stat.S_IRWXU)

home = os.environ['HOME']
cwd = os.getcwd()
status=False

shell = os.environ['SHELL'].split('/')[-1]
if shell == 'bash':
    if os.path.exists('%s/.bash_profile' % home):
        startupfile = '%s/.bash_profile'  % home
    elif os.path.exists('%s/.bash_login' % home):
        startupfile = '%s/.bash_login' % home
    elif os.path.exists('%s/.profile' % home):
        startupfile = '%s/.profile' % home
    else:
        startupfile = '%s/.bash_profile' % home
    status=True
elif shell == 'csh' or shell =='tcsh':
    startupfile = '%s/.cshrc' % home
    status=True
elif shell == 'zsh':
    startupfile = '%s/.zshrc' % home
    status=True
else:
    print('Can not find startup file to set envionment variables. \nPlease locate the startup file for your shell, and set manually\n For example for csh:\nsetenv PATH ${PATH}:/path/to_mirp/mirp/mirp\nsetenv PYTHONPATH /path/to_mirp/mirp/mirp')
    flag = False

if status:
    with open(startupfile, 'a') as f:
        f.write('\n\n### MiRP initialise ###\n')
        if shell == 'bash':
            f.write('export PATH=%s/mirp:$PATH\n' % cwd)
            f.write('export PYTHONPATH=%s/mirp:$PYTHONPATH\n' % cwd)
        elif shell == 'tcsh' or 'csh':
            f.write('setenv PATH ${PATH}:%s/mirp\n' % cwd)
            f.write('setenv PYTHONPATH %s/mirp\n' % cwd)
        elif shell == 'zsh':
            f.write('export PATH=%s/mirp:$PATH\n' % cwd)
            f.write('export PYTHONPATH=%s/mirp:$PYTHONPATH\n' % cwd)
    print('Install success! Please logout and in again for changes to take effect.')
else:
    print('Install failed.')

