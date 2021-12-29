import os

class figure_axis_xy(object):
        """define standart  XY Plot with reduced grafics"""

        def __init__(self,x_size=None,y_size=None,view_scale=None, size_tuple=None , fig_scale=None, container=False, dpi=180):
                import matplotlib.pyplot as plt
                import string
                import os

                if size_tuple is not None:
                    xsize, ysize =size_tuple[0], size_tuple[1]
                else:
                    xsize=x_size if x_size is not None else 8
                    ysize=y_size if y_size is not None else 5

                viewscale=view_scale if view_scale is not None else 0.5
                fig_scale=fig_scale if fig_scale is not None else 1

                self.label_letters=iter([i+') ' for i in list(string.ascii_lowercase)])

                if container:
                    self.fig=plt.figure(edgecolor='None',dpi=dpi*viewscale,figsize=(xsize*fig_scale, ysize*fig_scale),facecolor='w')
                else:
                    self.fig, self.ax=plt.subplots(num=None, figsize=(xsize*fig_scale, ysize*fig_scale), dpi=dpi*viewscale, facecolor='w', edgecolor='None')

        def make_clear_weak(self):
                #turn off axis spine to the right
                #self.fig.tight_layout()
                self.ax.spines['right'].set_color("none")
                self.ax.yaxis.tick_left() # only ticks on the left side
                self.ax.spines['top'].set_color("none")
                self.ax.xaxis.tick_bottom() # only ticks on the left side
        def make_clear(self):
                self.make_clear_weak()

        def make_clear_strong(self):
                #turn off axis spine to the right
                #self.fig.tight_layout()
                self.ax.spines['right'].set_color("none")
                self.ax.spines['left'].set_color("none")
                self.ax.yaxis.tick_left() # only ticks on the left side
                self.ax.spines['top'].set_color("none")
                self.ax.spines['bottom'].set_color("none")
                self.ax.xaxis.tick_bottom() # only ticks on the left side

        def tight(self):
                #turn off axis spine to the right
                self.fig.tight_layout()

        def label(self, x='x',y='y',t=None):

                self.ax.set_xlabel(x)
                self.ax.set_ylabel(y)
                self.ax.set_title(t, y=1.04)

        def save(self,name=None,path=None, verbose=True):
                import datetime
                import os

                savepath=path if path is not None else os.path.join(os.path.dirname(os.path.realpath('__file__')),'plot/')
                if not os.path.exists(savepath):
                    os.makedirs(savepath)
                #os.makedirs(savepath, exist_ok=True)
                name=name if name is not None else datetime.date.today().strftime("%Y%m%d_%I%M%p")
                #print(savepath)
                #print(name)
                extension='.pdf'
                full_name= (os.path.join(savepath,name)) + extension
                #print(full_name)
                self.fig.savefig(full_name, bbox_inches='tight', format='pdf', dpi=180)
                if verbose:
                    print('save at: '+name)

        def save_pup(self,name=None,path=None, verbose=True):
                import datetime
                import os
                import re

                name=re.sub("\.", '_', name)

                savepath=path if path is not None else os.path.join(os.path.dirname(os.path.realpath('__file__')),'plot/')
                if not os.path.exists(savepath):
                    os.makedirs(savepath)
                #os.makedirs(savepath, exist_ok=True)
                name=name if name is not None else datetime.date.today().strftime("%Y%m%d_%I%M%p")
                #print(savepath)
                #print(name)
                extension='.pdf'
                full_name= (os.path.join(savepath,name)) + extension
                #print(full_name)
                self.fig.savefig(full_name, bbox_inches='tight', format='pdf', dpi=300)
                if verbose:
                    print('save at: ',full_name)

        def save_light(self,name=None,path=None, verbose=True):
                import datetime
                import os

                savepath=path if path is not None else os.path.join(os.path.dirname(os.path.realpath('__file__')),'plot/')
                if not os.path.exists(savepath):
                    os.makedirs(savepath)
                #os.makedirs(savepath, exist_ok=True)
                name=name if name is not None else datetime.date.today().strftime("%Y%m%d_%I%M%p")
                #print(savepath)
                #print(name)
                extension='.png'
                full_name= (os.path.join(savepath,name)) + extension
                #print(full_name)
                self.fig.savefig(full_name, bbox_inches='tight', format='png', dpi=180)
                if verbose:
                    print('save with: ',name)


def write_log(hist, string, verbose=False, short=True , date=True):
    import datetime as datetime
    if short:
        now = datetime.datetime.now().strftime("%Y%m%d")
    else:
        now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    if date:
        message='\n'+now+' '+string
    else:
        message='\n '.ljust(16)+' '+string

    if verbose== True:
        print(message)
    elif verbose == 'all':
        print(hist+message)
    return hist+message

def json_load(name, path, verbose=False):
    import os
    import json
    full_name= (os.path.join(path,name+ '.json'))

    with open(full_name, 'r') as ifile:
        data=json.load(ifile)
    if verbose:
        print('loaded from: ',full_name)
    return data

def json_save(name, path, data, verbose=False, return_name=False):
    import json
    if not os.path.exists(path):
        os.makedirs(path)
    full_name_root=os.path.join(path,name)
    full_name= (os.path.join(full_name_root+ '.json'))
    with open(full_name, 'w') as outfile:
        json.dump(data, outfile)
    if verbose:
        print('save at: ',full_name)
    if return_name:
        return full_name_root
    else:
        return

def mkdirs_r(path):
    import os
    if not os.path.exists(path):
        os.makedirs(path)


def save_log_txt(name, path, hist,verbose=False):
    if not os.path.exists(path):
        os.makedirs(path)
    full_name= (os.path.join(path,name+ '.hist.txt'))
    with open(full_name, 'w') as ifile:
        ifile.write(str(hist))
    if verbose:
        print('saved at: ',full_name)
