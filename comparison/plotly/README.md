When a group of violin plots is specified, if there are many violin plots, each
ends up being too narrow to show the density clearly. To correct this, you can
change the `width` parameter. However, all the violin plots will be plotted
over each other, even when `violinmode = 'group'` rather than `overlay`. To
correct this, I hacked Plotly to permit the violins to be offset from each
other within a group.

To the URL of the page containing the plot, add `?intragroupoffset=0.3`, for
instance, to use this offset.

To use the patched Plotly, see https://github.com/plotly/plotly.js/blob/master/CONTRIBUTING.md#development.

In short:

    git clone https://github.com/plotly/plotly.js.git
    cd plotly.js
    npm install
    npm run pretest
    npm run bundle
    cp -a dist/plotly.min.js <dest>/
