import os
import sys

def get_env_info():
    env_path = os.path.dirname(sys.executable)
    env_name = os.path.basename(env_path)
    return env_name, env_path



def set_esmf_mk_path():
    env_name, env_path = get_env_info()
    if env_path:
        if env_path.endswith('/bin'): env_path= env_path[:-4]
        esmf_mk_path = os.path.join(env_path,'lib', 'esmf.mk')
        if os.path.isfile(esmf_mk_path):
            os.environ['ESMFMKFILE'] = esmf_mk_path
            print(f"Found esmf.mk at: {esmf_mk_path}")
        else:
            print('esmf.mk not found in the specified conda environment.')
    else:
        print('Unable to determine the Anaconda environment.')
        
set_esmf_mk_path()
    
        
        
