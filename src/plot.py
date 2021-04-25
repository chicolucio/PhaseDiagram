from matplotlib import pyplot as plt
from phase_diagram import ureg


class Plot:
    def __init__(self, x_unit, y_unit, x_label='', y_label='', ax=None, scale_log=True, legend=False, title=True,
                 title_text=''):
        """
        Parameters
        ----------
        x_unit : str
            pint unit of the independent variable
        y_unit : str
            pint unit of the dependent variable
        x_label : str
            label of the independent variable
        y_label : str
            label of the dependent variable
        ax : matplotlib axis, optional
            axis where the plot will be shown. If None, one will be created
        scale_log : bool, default=True
            if the y-axis will have a log scale
        legend : bool, default=False
            if a legend will be shown
        title : bool, default=True
            if the plot will have a title
        title_text : str, default=''
            title text
        """
        self.x_unit = x_unit
        self.y_unit = y_unit
        self.x_label = x_label
        self.y_label = y_label
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

        if self.scale_log:
            self.ax.set_yscale('log')
            self.ax.set_ylabel('log({} / {:~P})'.format(self.y_label, ureg(self.y_unit).units))
        else:
            # setting the y-axis to scientific notation and
            # getting the order of magnitude
            self.ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
            self.ax.yaxis.major.formatter._useMathText = True
            self.ax.figure.canvas.draw()  # Update the text
            order_magnitude = self.ax.yaxis.get_offset_text().get_text().replace('\\times', '')
            self.ax.yaxis.offsetText.set_visible(False)

            self.ax.set_ylabel('{} / '.format(self.y_label) + order_magnitude + ' {:~P}'.format(ureg(self.y_unit).units))

        self.ax.set_xlabel('{} / {:~P}'.format(self.x_label, ureg(self.x_unit).units))

        if self.legend:
            self.ax.legend(loc='best', fontsize=14)

        if self.title:
            self.ax.set_title(self.title_text, fontsize=18)

        return self.ax

    def plot_arrays(self, tuple_two_arrays, limit=None, label='', **kwargs):
        """
        Creates a plot based on two arrays

        Parameters
        ----------
        tuple_two_arrays : tuple
            tuple with two arrays (x array, y array) with pint units
        limit : float
            creates a boolean mask in y array, limit y to y< limit
        label : str
            label in the legend
        **kwargs : optional
            matplotlib arguments
        """
        x, y = tuple_two_arrays
        try:
            y = y[y < limit]
            x = x[:len(y)]
        except:
            pass
        self.ax.plot(x.to(self.x_unit), y.to(self.y_unit), label=label, **kwargs)
        self.plot_customization()

    def plot_point(self, tuple_point, label='', **kwargs):
        """
        Creates a point in a plot based

        Parameters
        ----------
        tuple_point : tuple
            tuple with two values (x value, y value) with pint units
        label : str
            label in the legend
        **kwargs : optional
            matplotlib arguments
        """
        self.ax.scatter(tuple_point[0].to(self.x_unit), tuple_point[1].to(self.y_unit), label=label, **kwargs)
        self.plot_customization()
