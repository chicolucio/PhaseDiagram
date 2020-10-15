import matplotlib.pyplot as plt
from . import ureg


def _plot_params(ax=None):
    """Internal function for plot parameters.
    Parameters
    ----------
    ax : Matplotlib axes, optional
        axes where the graph will be plotted, by default None
    """
    linewidth = 2
    size = 12

    # grid and ticks settings
    ax.minorticks_on()
    ax.grid(b=True, which='major', linestyle='--',
            linewidth=linewidth - 0.5)
    ax.grid(b=True, which='minor', axis='both',
            linestyle=':', linewidth=linewidth - 1)
    ax.tick_params(which='both', labelsize=size + 2)
    ax.tick_params(which='major', length=6, axis='both')
    ax.tick_params(which='minor', length=3, axis='both')

    # labels and size
    ax.xaxis.label.set_size(size + 4)
    ax.yaxis.label.set_size(size + 4)
    # ax.title.set_fontsize(size+6)  # not working, don't know why...

    return


def plot_arrays(tuple_two_arrays, tuple_units,  limit=None, ax=None, **kwargs):
    temperature, pressure = tuple_two_arrays
    try:
        pressure = pressure[pressure < limit]
        temperature = temperature[:len(pressure)]
    except TypeError:
        pass
    ax.plot(temperature.to(tuple_units[0]), pressure.to(tuple_units[1]), **kwargs)


def plot_point(tuple_point, ax=None, **kwargs):
    ax.scatter(tuple_point[0], tuple_point[1], **kwargs)


def plot_customization(ax=None, T_unit='K',
         P_unit='Pa', scale_log=True, legend=False, title=True,
         title_text=''):

    if scale_log:
        ax.set_yscale('log')
        ax.set_ylabel('log(Pressure / {:~P})'.format(ureg(P_unit).units))
    else:
        # setting the y-axis to scientific notation and
        # getting the order of magnitude
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        ax.yaxis.major.formatter._useMathText = True
        ax.figure.canvas.draw()  # Update the text
        order_magnitude = ax.yaxis.get_offset_text().get_text().replace('\\times', '')
        ax.yaxis.offsetText.set_visible(False)

        ax.set_ylabel('Pressure / ' + order_magnitude +
                      ' {:~P}'.format(ureg(P_unit).units))

    ax.set_xlabel('Temperature / {:~P}'.format(ureg(T_unit).units))

    if legend:
        ax.legend(loc='best', fontsize=14,
                  title=self.format_formula(), title_fontsize=14)

    if not title:
        pass
    elif title_text == '':
        ax.set_title('Calculated phase diagram - ' + self.format_formula(),
                     fontsize=18)
    else:
        ax.set_title(title_text, fontsize=18)

    return ax
