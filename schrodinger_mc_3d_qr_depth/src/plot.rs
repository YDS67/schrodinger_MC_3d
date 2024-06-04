use plotly::color::{NamedColor, Rgb};
use plotly::common::{Anchor, Font, Line, Marker, MarkerSymbol, Mode, Title};
use plotly::layout::{Axis, Legend, Shape, ShapeLine, ShapeType, ItemSizing, Margin};
use plotly::{ImageFormat, Layout, Plot, Scatter};

pub struct PlotPar{
    pub xlab: String,
    pub ylab: String,
    pub title: String,
    pub flnm: String,
    pub legends: Vec<String>
}

impl PlotPar{
    pub fn new(xlab: &str, ylab: &str, title: &str, flnm: &str, legends: Vec<String>) -> PlotPar {
        PlotPar {
            xlab: format!("{}", xlab),
            ylab: format!("{}", ylab),
            title: format!("{}", title),
            flnm: format!("{}", flnm),
            legends,
        }
    }
}

pub fn line_and_scatter_plot(x1: &Vec<f64>, y1: &Vec<f64>, x2: &Vec<f64>, y2: &Vec<f64>, plot_par: &PlotPar) {
    let bgcol = Rgb::new(255, 255, 255);
    let linecol1 = NamedColor::DarkBlue;
    let linecol2 = NamedColor::DarkRed;
    let forecol = Rgb::new(0, 0, 0);
    let gridcol = Rgb::new(180, 180, 180);
    let transp = NamedColor::Transparent;
    let thick: usize = 3;
    let medium: usize = 3;
    let _thin: usize = 2;
    let msize: usize = 10;
    let fsz_title: usize = 35;
    let fsz_legend: usize = 35;
    let fsz_ticks: usize = 30;
    let fsz_axes: usize = 35;

    let trace1 = Scatter::new(x1.clone(), y1.clone())
        .name(&plot_par.legends[0])
        .mode(Mode::Lines)
        .line(Line::new().color(linecol1).width(medium as f64)
        //.marker(Marker::new().size(msize).symbol(MarkerSymbol::Circle),
    );
        
    let trace2 = Scatter::new(x2.clone(), y2.clone())
        .name(&plot_par.legends[1])
        .mode(Mode::Markers)
        //.line(Line::new().color(linecol2).width(medium as f64))
        .marker(Marker::new().size(msize).color(linecol2).symbol(MarkerSymbol::Circle)
    );

    let title = Title::new(&plot_par.title)
        .font(Font::new().size(fsz_title).family("Serif").color(forecol));

    let legend = Legend::new()
        .x(0.99)
        .x_anchor(Anchor::Right)
        .y(1.0 - 0.0133)
        .y_anchor(Anchor::Top)
        .font(Font::new().size(fsz_legend).color(forecol).family("Serif"))
        .border_width(medium)
        .border_color(forecol)
        .background_color(bgcol)
        .item_width(52)
        .item_sizing(ItemSizing::Trace);

    let axis = Axis::new()
        .position(0.0)
        .show_line(true)
        .line_color(forecol)
        .line_width(thick)
        .tick_length(9)
        .tick_width(medium)
        .tick_color(forecol)
        .tick_font(Font::new().color(forecol))
        .zero_line(false)
        .show_grid(true)
        .grid_color(gridcol);

    let axisx = axis.clone().title(
        Title::new(&plot_par.xlab)
            .font(Font::new().size(fsz_axes).color(forecol).family("Serif")));

    let axisy = axis
        .clone()
        .title(Title::new(&plot_par.ylab)
            .font(Font::new().size(fsz_axes).color(forecol).family("Serif")))
        .tick_angle(270.0);

    let line_top = Shape::new()
        .shape_type(ShapeType::Line)
        .x_ref("paper")
        .y_ref("paper")
        .x0(0.)
        .y0(1.)
        .x1(1.)
        .y1(1.)
        .line(ShapeLine::new().color(forecol).width(thick as f64));

    let line_right = Shape::new()
        .shape_type(ShapeType::Line)
        .x_ref("paper")
        .y_ref("paper")
        .x0(1.)
        .y0(0.)
        .x1(1.)
        .y1(1.)
        .line(ShapeLine::new().color(forecol).width(thick as f64));

    let mut layout = Layout::new()
        .width(1024)
        .height(768)
        .font(Font::new().size(fsz_ticks))
        .title(title)
        .legend(legend)
        .show_legend(true)
        .x_axis(axisx)
        .y_axis(axisy)
        .plot_background_color(transp)
        .paper_background_color(bgcol)
        .margin(Margin::new().left(105).bottom(105));

    layout.add_shape(line_top);
    layout.add_shape(line_right);

    let mut plot = Plot::new();
    plot.add_trace(trace1);
    plot.add_trace(trace2);
    plot.set_layout(layout);

    //plot.write_html(flnm);
    plot.write_image(&plot_par.flnm, ImageFormat::PDF, 1280, 960, 1.0);
    //plot.write_image(flnm, ImageFormat::PNG, 1280, 960, 1.0);
}
