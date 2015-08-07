# coding: utf-8

import csv
import numpy as np
import os

import matplotlib.pyplot as plt
import matplotlib
import operator

from scipy.stats import linregress

def plot_original_fc_over_length():
    gene_fc = []
    gene_length = []
    gene_id = []

    #are we reading the  first line from the file
    first_row = True

    with open('Test_IBA.txt', 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        for row in reader:

            #lets skip the first row, it is just the headers
            if first_row:
                first_row = False
                continue

            id, fc, length = row

            if 'N/A' not in length:
                try:
                    #cast the numbers from strings to floats
                    fc_float = float(fc)
                    length_float = float(length)

                    gene_id.append(id)
                    gene_fc.append(fc_float)
                    gene_length.append(length_float)

                except Exception as e:
                    pass

    #convert the lists to arrays, lots of mathematical operations
    #much easier to apply to arrays
    gene_length = np.array(gene_length)
    gene_fc = np.array(gene_fc)

    log_gene_length = np.log10(gene_length)

    plt.scatter(log_gene_length, gene_fc.clip(-10, 10))
    plt.show()
    import IPython
    IPython.embed()
    assert False


def runningMean(fc, length, step_size, window_size, use_50_percentile=False):
    out = np.array([(y, x) for (y, x) in sorted(zip(length, fc))])
    fc_sorted_by_length = out[:, 1]
    sorted_length = out[:, 0]

    x = []
    y = []

    for ctr in range(0, len(fc_sorted_by_length)-window_size, step_size):

        fcs = fc_sorted_by_length[ctr:(ctr+window_size)]
        lengths = sorted_length[ctr:ctr+window_size]

        #sort by fc and remove the bottom25% and top 25%
        #in order to smooth the plots
        if use_50_percentile:
            fcs_lengths = zip(fcs, lengths)
            fcs_lengths.sort(key=operator.itemgetter(0))
            fcs_lengths = np.array(fcs_lengths)

            fcs = fcs_lengths[window_size/4:(window_size-window_size/4), 0] * 2.0
            lengths = fcs_lengths[window_size/4:(window_size-window_size/4), 1] * 2.0

        y.append(np.sum(fcs))
        x.append(np.sum(lengths))

    x = np.array(x)
    y = np.array(y)

    return x/float(window_size), y/float(window_size)


def get_fc_and_length_from_file(filename):

    gene_fc = []
    gene_length = []
    gene_id = []

    #are we reading the  first line from the file
    first_row = True
    #na_count = 0

    with open(filename, 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        for row in reader:

            #lets skip the first row, it is just the headers
            if first_row:
                first_row = False
                continue

            id, fc, length = row

            if fc == '' and length == '' and id == '':
                continue

            if "N/A" in length:
                continue

            gene_id.append(id)

            if float(fc) < 0:
                fc = -1/float(fc)

            gene_fc.append(float(fc))
            gene_length.append(float(length))

    gene_length = np.array(gene_length)
    gene_fc = np.array(gene_fc)

    return gene_id, gene_fc, gene_length


def calculate_fc_from_file(filename):

    value_1 = []
    value_2 = []
    gene_length = []
    gene_id = []

    #are we reading the  first line from the file
    first_row = True
    #na_count = 0

    with open(filename, 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        for row in reader:

            #lets skip the first row, it is just the headers
            if first_row:
                first_row = False
                continue

            id, v1, v2, length = row

            if v1 == '' and length == '' and id == '':
                continue

            if "N/A" in length:
                continue

            # if float(v1) < 0.0000001:
            #     #print v1, v2
            #     continue
            #
            # if float(v2) < 0.0000001:
            #     #print v1, v2
            #     continue

            gene_id.append(id)
            value_1.append(float(v1))
            value_2.append(float(v2))
            gene_length.append(float(length))

    gene_length = np.array(gene_length)
    value_1 = np.array(value_1)
    value_2 = np.array(value_2)

    gene_fc = (value_1 + .00001) / (value_2 + .00001)
    #gene_fc = (value_1 ) / (value_2 )

    return gene_id, gene_fc, gene_length


def run_file(filename="data/12hr_test.txt",
             output_directory="plots/",
             step_size=40,
             window_size=200,
             save=True,
             use_50_percentile=False):

    #if the file does not have an fc column, then we need to calculate it
    try:
        gene_id, gene_fc, gene_length = calculate_fc_from_file(filename)
    #if the file has the fc column, we can just load it
    except ValueError as e:
        gene_id, gene_fc, gene_length = get_fc_and_length_from_file(filename)

    log_gene_length = np.log10(gene_length)
    log_gene_fc = np.log2(gene_fc)

    running_mean_log_lengths, running_mean_log_fc = runningMean(log_gene_fc,
                                                                log_gene_length,
                                                                step_size,
                                                                window_size,
                                                                use_50_percentile=use_50_percentile)

    running_mean_lengths_in_kb = np.power(10, running_mean_log_lengths)/1000.0
    running_mean_log_lengths_in_kb = running_mean_log_lengths - 3.0

    plot_regression = step_size is 1 and win_size is 1
    #
    if plot_regression:
        slope, intercept, r_value, p_value, std_err = linregress(running_mean_lengths_in_kb, running_mean_log_fc)

        print r_value
        print p_value
        print "r_value: %f\np_value: %f\nstd_err: %f\n" % (r_value, p_value, std_err)

    #clear the old figure before adding a new plot
    plt.clf()

    ##########################################
    #Plot parameters
    ##########################################
    # width of the entire figure in inches
    WIDTH_INCHES = 2*4
    # height of the entire figure in inches.
    HEIGHT_INCHES = 1.5*4
    TITLE_FONTSIZE = '28'
    # font size for 0,10,20, etc on each axis
    TICK_LABEL_FONT_SIZE = 14
    # Gene Length font size
    AXES_LABEL_FONT_SIZE = '28'
    # dots per inch.
    DPI = 4*72
    # how big are the little circles that are plotted on the graph.
    PLOT_MARKER_SIZE = 20
    # if the ylabel is getting cut off, make this larger.
    LEFT_ADJUST = 0.18
    # if the x label is getting cut off, make this larger.
    BOTTOM_ADJUST = 0.18
    ############################################

    fig = plt.figure(figsize=(WIDTH_INCHES, HEIGHT_INCHES))

    #set the plot to be log scale
    if not plot_regression:
        plt.xscale('log', basex=10)

    scatter = plt.scatter(running_mean_lengths_in_kb, running_mean_log_fc, s=PLOT_MARKER_SIZE)
    #scatter = plt.scatter(running_mean_log_lengths_in_kb, running_mean_log_fc)

    #add a horizontal line at y = 0, 'k' is the variable for black colored line
    plt.axhline(y=0, color='k')
    #
    if plot_regression:
         plt.plot(running_mean_lengths_in_kb, slope*running_mean_lengths_in_kb + intercept, color="r")
    #
    #generate the title for the figure from the name of the file, step_size, and window_size
    fig_title = filename.split('/')[-1].replace(".txt", " ").replace("_", " ") + str(step_size) + " " + str(window_size)
    if use_50_percentile:
        fig_title += " 50 Percentile"

    if plot_regression:
         fig_title += " Regression"

    from matplotlib import font_manager as fm

    title_font = {'size': TITLE_FONTSIZE, 'color': 'black', 'weight': 'normal', 'verticalalignment': 'bottom'}
    prop = fm.FontProperties(fname='/usr/share/fonts/truetype/msttcorefonts/Arial.ttf')

    for tick in scatter.axes.get_xmajorticklabels():
        tick.set_fontsize(TICK_LABEL_FONT_SIZE)

    for tick in scatter.axes.get_ymajorticklabels():
        tick.set_fontsize(TICK_LABEL_FONT_SIZE)

    scatter.axes.set_xlabel("Gene Length (kb)",  fontsize=AXES_LABEL_FONT_SIZE)
    scatter.axes.set_ylabel('Mean Log$_2$ Fold Change', fontsize=AXES_LABEL_FONT_SIZE)

    scatter.axes.tick_params(axis='both', which='major', pad=15)
    scatter.axes.tick_params(axis='y', direction='out')
    scatter.axes.tick_params(axis='x', direction='out')

    scatter.axes.get_xaxis().tick_bottom()#removes ticks from top
    scatter.axes.get_yaxis().tick_left()#removes ticks on rhs
    scatter.axes.spines["right"].set_visible(False)
    scatter.axes.spines["top"].set_visible(False)

    scatter.axes.set_title(fig_title, fontproperties=prop, **title_font)

    plt.gcf().subplots_adjust(bottom=BOTTOM_ADJUST, left=LEFT_ADJUST)

    #want the x label to be 1,10,100, ... not 10^1, 10^2,...
    scatter.axes.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

    if save:
        plt.savefig(output_directory + fig_title.replace(" ", "_") + ".pdf", dpi=DPI)
    else:
        plt.show()

######
#THIS IS THE START OF OUR PROGRAM
######
if __name__ == "__main__":
    #
    # # where we are getting data from
    # DATA_DIRECTORY = "work/"
    # #where we are saving the plots
    # OUTPUT_DIRECTORY = "work_plots/"
    # #whether to save the plots, if false, they are displayed instead.
    # SAVE = False
    # #we have the fold change in the .txt file, don't need to calculate it
    # CALCULATE_FC = False
    #
    # SLIDING_WINDOW_STEP_SIZE = 1
    # SLIDING_WINDOW_SIZE = 80

    # where we are getting data from
    DATA_DIRECTORY = "data/"
    #where we are saving the plots
    OUTPUT_DIRECTORY = "plots/"
    #whether to save the plots, if false, they are displayed instead.
    SAVE = True
    SLIDING_WINDOW_STEP_SIZE = 10
    SLIDING_WINDOW_SIZE = 800

    #for each sliding window, average the
    #whole window, or only use middle 50 percentile
    USE_50_PERCENTILE = False

    if not os.path.exists(OUTPUT_DIRECTORY):
        os.mkdir(OUTPUT_DIRECTORY)

    # #USE THIS FOR REGRESSION
    win_sizes = [1]
    step_sizes = [1]

    for filename in os.listdir(DATA_DIRECTORY):
        for win_size in win_sizes:
            for step_size in step_sizes:
                print
                print "working on: " + str(filename)
                print "win size: " + str(win_size)
                print "step size: " + str(step_size)

                run_file(filename=DATA_DIRECTORY + filename,
                         output_directory=OUTPUT_DIRECTORY,
                         step_size=step_size,
                         window_size=win_size,
                         save=SAVE,
                         use_50_percentile=USE_50_PERCENTILE)