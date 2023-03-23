import numpy as np
import os

class Model(object):
    def __init__(self, filename, file_format = None):
        if filename is None or not os.path.isfile(filename):
            raise Exception("Given filename was not found!")

        if file_format is None and filename[-3:] == '.nd':
            file_format = 'nd'

        if file_format is None and filename[-5:] == '.tvel':
            file_format = 'tvel'
        
        if file_format is None:
            raise Exception('Cannot find file format')

        self.file_format = file_format
        self.filename = filename
        self.modified = False

        self.__model__, self.__is_flat, self.__has_ocean__ = Model.load_nd(filename) if self.file_format == 'nd' else Model.load_tvel(filename)

    @staticmethod
    def __new_model__():
        
        return {
            'z' : [],
            'vp' : [],
            'vs' : [],
            'rho' : [],
            'qal' : [],
            'qemu' : []
            }

    @staticmethod
    def load_nd(filename):
        model = Model.__new_model__()
        with open(filename, "r") as fio:
            for line in fio:
                if line.strip() in ['mantle', 'outer-core', 'inner-core']: continue
                z, vp, vs, dens, qal, qemu = line.split()
                model['z'].append(float(z))
                model['vp'].append(float(vp))
                model['vs'].append(float(vs))
                model['rho'].append(float(dens))
                model['qal'].append(float(qal))
                model['qemu'].append(float(qemu))

        model['z']    = np.array(model['z'])
        model['vp']   = np.array(model['vp'])
        model['vs']   = np.array(model['vs'])
        model['rho']  = np.array(model['rho'])
        model['qal']  = np.array(model['qal'])
        model['qemu'] = np.array(model['qemu'])

        return (model, False, (True if model['vs'][0] == 0.0 else False))

    @staticmethod
    def load_tvel(filename):
        model = Model.__new_model__()
        
        with open(filename, "r") as fio:
            for i,line in enumerate(fio):
                if i < 2: continue
                z, vp, vs, dens = line.split()
                model['z'].append(float(z))
                model['vp'].append(float(vp))
                model['vs'].append(float(vs))
                model['rho'].append(float(dens))
        
        model['z']   = np.array(model['z'])
        model['vp']  = np.array(model['vp'])
        model['vs']  = np.array(model['vs'])
        model['rho'] = np.array(model['rho'])

        model['qal']  = None
        model['qemu'] = None
        
        return (model, False, (True if model['vs'][0] == 0.0 else False))

    def __str__(self):
        txt = "z\tvp\tvs\trho\tqal\tqemu\n"
        try:
            for z,vp,vs,rho,qal,qemu in zip(*self.__model__.values()):
                txt += f'{z:.1f}\t{vp:.4f}\t{vs:.4f}\t{rho:.2f}\t{qal}\t{qemu}\n'
        except TypeError:
            for z,vp,vs,rho in zip(*list(self.__model__.values())[:4]):
                txt += f'{z:.1f}\t{vp:.4f}\t{vs:.4f}\t{rho:.2f}\t-\t-\n'            
        return txt

    def get(self, var):
        if not (var in self.__model__) or (self.__model__[var] is None):
            raise KeyError('Variable not in model.')

        return self.__model__[var].copy()

    def thick_get(self, variables, max_depth = None,
                  handle_water = True, get_z_values = False):
        
        if variables is None: KeyError('Bad variable to get.')
        if not isinstance(variables, list): variables = [ variables ]
        if 'z' in variables: raise KeyError('Z is not allowed.')
        
        for onevar in variables:
            if onevar not in [ 'vp', 'vs', 'rho' ]:
                raise Exception(f"Bad variable name = {onevar}")
        
        items = []

        z     = self.get('z')
        thick = z[1:]-z[:-1]

        # Remove the 0 thickness layers
        ff = (thick != 0)
        if max_depth is not None:
            ff = ff&(z[:-1] <= max_depth)
        
        # Collect needed variables
        for varn in variables:
            var = self.get(varn)
            
            var = (var[1:]+var[:-1])/2
            var = var[ff]
            
            # Check vs is zero somewhere
            if handle_water and varn == 'vs':
                var[var < 0.001] = 1E-5
            
            # Save variable
            items.append(var)

        items.insert(0,thick[ff])

        # Usefull for check plot
        if get_z_values:
            variables.append('depth')
            items.append(z[:-1][ff])
        
        return tuple(items)

    def insert_prem_ocean(self, depth = 3.0):
        if self.__has_ocean__:
            return self.__has_ocean__

        if depth >= self.__model__['z'][1]:
            raise Exception("Cannot insert an ocean deeper than first layer")

        # Fix top
        self.__model__['z'][0] = depth
        self.__model__['z']   = np.insert(self.__model__['z']  , 0, [0.0, depth])
        self.__model__['vp']  = np.insert(self.__model__['vp'] , 0, [1.45, 1.45])
        self.__model__['vs']  = np.insert(self.__model__['vs'] , 0, [0.0, 0.0])
        self.__model__['rho'] = np.insert(self.__model__['rho'], 0, [1.02, 1.02])

        if self.__model__['qal'] is not None:
            self.__model__['qal'] = np.insert(self.__model__['qal'], [0,1], [57822, 57822])

        if self.__model__['qemu'] is not None:
            self.__model__['qemu'] = np.insert(self.__model__['qemu'], [0,1], [0, 0])

        self.__has_ocean__ = True

        return self.__has_ocean__

    def save_model96(self, filename, cmt = "None"):
        sph_type = "FLAT EARTH"

        if not self.__is_flat:
            sph_type = "SPHERICAL EARTH"

        heads = [ "H", "VP", "VS", "RHO", "QP", "QS", "ETAP", "ETAS", "FREFP", "FREFS" ]

        with open(filename, "w") as fio:
            print("MODEL.01", file = fio)
            print(cmt, file = fio)
            print("ISOTROPIC", file = fio)
            print("KGS", file = fio)
            print(sph_type, file = fio)
            print("1-D", file = fio)
            print("CONSTANT VELOCITY", file = fio)
            print("LINE08", file = fio)
            print("LINE09", file = fio)
            print("LINE10", file = fio)
            print("LINE11", file = fio)
            print("\t".join(heads), file = fio)
            h,vp,vs,rho = self.thick_get(["vp", "vs", "rho"])
            for ih, ivp, ivs, irho in zip(h, vp, vs, rho):
                print(f'{ih:.1f}\t{ivp:.4f}\t{ivs:.4f}\t{irho:.2f}\t0.0\t0.0\t0.0\t0.0\t1.0\t1.0', file = fio)

        return True

    def layer_mask(self, what):
        if self.modified: raise Exception("Cannot create mask on a modified model")
        
        z = self.get('z')
        i = np.arange(len(z))
        
        if what == 'c':
            mask = (z<=24.4)
            mask[i[mask][-1]] = False

        elif what == 'l':
            mask = (z>=24.4)&(z<=220)
            mask[i[mask][[0,-1]]] = False
        
        elif what == 'u':
            mask = (z>=220.)&(z<=400)
            mask[i[mask][[0,-1]]] = False

        elif what == 'tz':
            mask = (z>=400.)&(z<=670)
            mask[i[mask][[0,-1]]] = False

        elif what == 'lm':
            mask = (z>=670.)&(z<=2891)
            mask[i[mask][[0,-1]]] = False

        elif what == 'oc':
            mask = (z>=2891.)&(z<=5149.5)
            mask[i[mask][[0,-1]]] = False

        elif what == 'ic':
            mask = (z>=5149.5)
            mask[i[mask][0]] = False

        else:
            raise Exception(f"Bad layer {what}")

        return mask

    @staticmethod
    def T(z, v, ll, lu, scale_vel):
        if ll is None:
            ll = z[0]
        
        if lu is None:
            lu = z[-1]
        
        if ll >= z[-1]: raise Exception("Wow! ll is too deep!")
        if lu <= z[0]: raise Exception("Wow! lu is too shallow!")
        
        vl = np.interp(ll, z, v)
        vu = np.interp(lu, z, v)
        
        if scale_vel:
            if ll < z[0]:
                vl = v[0]-(z[0]-ll)*((v[1]-v[0])/(z[1]-z[0]))

            if lu > z[-1]:
                vu = v[-1]-(z[-1]-lu)*((v[-1]-v[-2])/(z[-1]-z[-2]))
        
        # Depth Change
        ##
        z = z.copy()
        z -= z.min()
        z /= z.max()
        z *= (lu-ll)
        z += ll

        # Vel Change
        ##
        v = v.copy()
        if scale_vel:
            v -= v[0]
            v /= v[-1]
            v *= (vu-vl)
            v += vl
        
        return z,v

    def mod(self, controls, inplace = False):
        if self.modified:
            raise Exception("Cannot modify a modified model")
        
        if not inplace:
            _m = Model(self.filename, self.file_format)
            
            if self.__has_ocean__:
                _m.insert_prem_ocean(self.__model__['z'][1])
            
            _m.mod(controls, True)
            
            return _m
        
        z = np.array([])
        vp = np.array([])
        vs = np.array([])

        for l in ['c', 'l', 'u', 'tz', 'lm', 'oc', 'ic']:
            zl   = self.get('z')[self.layer_mask(l)]
            vpl  = self.get('vp')[self.layer_mask(l)]
            vsl  = self.get('vs')[self.layer_mask(l)]
            
            if l in controls:
                _, vpl = Model.T(zl, vpl, *controls[l])
                zl, vsl = Model.T(zl, vsl, *controls[l])
            
            z  = np.append(z, zl)
            vp = np.append(vp, vpl)
            vs = np.append(vs, vsl)
        
        # Well, some discontinuities dissapear on rho
        cur_z = self.get('z')
        cur_rho = self.get('rho')
        rho = np.interp(z, cur_z, cur_rho)
        mask = z != cur_z
        
        # Load up
        self.__model__['z']   = z
        self.__model__['vp']  = vp
        self.__model__['vs']  = vs
        self.__model__['rho'][mask] = rho[mask]
        self.__model__['qal'] = None
        self.__model__['qemu'] = None
        
        return None


if __name__ == "__main__":
    import sys
    
    prem = Model("/usr/local/TauP-2.1.1/StdModels/prem.nd")
    
    
    # ~ prem.save_model96("prem.model96")
    # ~ prem.insert_prem_ocean()
    # ~ prem.save_model96("prem-ocean.model96")
    
    sys.exit(0)
