from matplotlib import pyplot as plt

from phase_diagram import ureg


class Plot:
    def __init__(self, x_unit, y_unit, ax=None, scale_log=True, legend=False, title=True, title_text=''):
        self.x_unit = x_unit
        self.y_unit = y_unit
        self.ax = ax
        self.scale_log = scale_log
        self.legend = legend
        self.title = title
        self.title_text = title_text

        if self.ax is None:
            fig, self.ax = plt.subplots(figsize=(10, 8), facecolor=(1.0, 1.0, 1.0))

    def plot_customization(self):

        linewidth = 2
        size = 12

        # grid and ticks settings
        self.ax.minorticks_on()
        self.ax.grid(b=True, which='major', linestyle='--',
                     linewidth=linewidth - 0.5)
        self.ax.grid(b=True, which='minor', axis='both',
                     linestyle=':', linewidth=linewidth - 1)
        self.ax.tick_params(which='both', labelsize=size + 2)
        self.ax.tick_params(which='major', length=6, axis='both')
        self.ax.tick_params(which='minor', length=3, axis='both')

        # labels and size
        self.ax.xaxis.label.set_size(size + 4)
        self.ax.yaxis.label.set_size(size + 4)
        # ax.title.set_fontsize(size+6)  # not working, don't know why...

        if self.scale_log:
            self.ax.set_yscale('log')
            self.ax.set_ylabel('log(Pressure / {:~P})'.format(ureg(self.y_unit).units))
        else:
            # setting the y-axis to scientific notation and
            # getting the order of magnitude
            self.ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
            self.ax.yaxis.major.formatter._useMathText = True
            self.ax.figure.canvas.draw()  # Update the text
            order_magnitude = self.ax.yaxis.get_offset_text().get_text().replace('\\times', '')
            self.ax.yaxis.offsetText.set_visible(False)

            self.ax.set_ylabel('Pressure / ' + order_magnitude +
                          ' {:~P}'.format(ureg(self.y_unit).units))

        self.ax.set_xlabel('Temperature / {:~P}'.format(ureg(self.x_unit).units))

        if self.legend:
            self.ax.legend(loc='best', fontsize=14,
                           title=self.format_formula(), title_fontsize=14)

        if not self.title:
            pass
        elif self.title_text == '':
            # TODO: resolver o format_formula
            # self.ax.set_title('Calculated phase diagram - ' + self.format_formula(),
            #                   fontsize=18)
            pass
        else:
            self.ax.set_title(self.title_text, fontsize=18)

        return self.ax

    def plot_arrays(self, tuple_two_arrays, limit=None, **kwargs):
        x, y = tuple_two_arrays
        try:
            y = y[y < limit]
            x = x[:len(y)]
        except:
            pass
        self.ax.plot(x.to(self.x_unit), y.to(self.y_unit), **kwargs)
        self.plot_customization()

    def plot_point(self, tuple_point, **kwargs):
        self.ax.scatter(tuple_point[0].to(self.x_unit), tuple_point[1].to(self.y_unit), **kwargs)
        self.plot_customization()