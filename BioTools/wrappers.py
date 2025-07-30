import matplotlib.pyplot as plt
from functools import wraps
from os import path, makedirs

def savefig(default_filename, default_rel_path=None, **default_savefig_kwargs):
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            filename = kwargs.pop('filename', default_filename)
            rel_path = kwargs.pop('rel_path', default_rel_path)
            format = kwargs.get('format', 'png')
            save = kwargs.pop('save', False)
            savefig_kwargs = default_savefig_kwargs.copy()
            savefig_kwargs.update(kwargs.pop('savefig_kwargs', {}))

            savefig_kwargs.setdefault('bbox_inches', 'tight')
            savefig_kwargs.setdefault('dpi', 300)
            savefig_kwargs.setdefault('format', 'png')
            
            
            result = func(*args, **kwargs)
            # Получаем расширение файла
            _, ext = path.splitext(filename)
            ext = ext.lower().lstrip('.')
            # Если расширение указано явно, используем его как формат
            if ext:
                savefig_kwargs['format'] = ext
                save_filename = filename
            else:
                save_filename = f"{filename}.{format}"
                print(f'{save_filename=}')

            if savefig_kwargs['format'] == 'tiff':
                savefig_kwargs.setdefault('pil_kwargs', {"compression": "tiff_lzw"})


            if save:
                save_path = save_filename
                print(f'{save_path=}')
                if rel_path:
                    if not path.exists(rel_path):
                        makedirs(rel_path)
                    save_path = path.join(rel_path, save_filename)
                    print(f'{save_path=}')
                plt.savefig(save_path, **savefig_kwargs)
                print(f"Plot saved to file: {save_path}")
            plt.show()
            return result
        return wrapper
    return decorator