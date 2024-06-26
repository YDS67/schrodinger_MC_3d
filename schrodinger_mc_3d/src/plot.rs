use plotly::color::{NamedColor, Rgb};
use plotly::common::{Anchor, Font, Line, Marker, MarkerSymbol, Mode, Title};
use plotly::layout::{Axis, Legend, Shape, ShapeLine, ShapeType, ItemSizing, Margin};
use plotly::{ImageFormat, Layout, Plot, Scatter};

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum LegendAl{
    BottomRight,
    TopRight,
    BottomLeft,
    TopLeft,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum LineOrPoints{
    Line,
    Points,
    LineAndPoints,
}

pub struct PlotPar{
    pub xlab: String,
    pub ylab: String,
    pub title: String,
    pub flnm: String,
    pub legends: Vec<String>,
    pub legend_al: LegendAl,
    pub line_or_points: Vec<LineOrPoints>,
}

impl PlotPar{
    pub fn new(xlab: &str, ylab: &str, title: &str, flnm: &str, legends: Vec<String>) -> PlotPar {
        PlotPar {
            xlab: format!("{}", xlab),
            ylab: format!("{}", ylab),
            title: format!("{}", title),
            flnm: format!("{}", flnm),
            legends,
            legend_al: LegendAl::TopRight,
            line_or_points: vec![LineOrPoints::Line; 20],
        }
    }
}

pub fn line_plot(x: &Vec<Vec<f64>>, y: &Vec<Vec<f64>>, plot_par: &PlotPar) {
    let lines_number = x.len();
    let bgcol = Rgb::new(255, 255, 255);
    let cols = vec![
        NamedColor::MidnightBlue, 
        NamedColor::Maroon,
        NamedColor::Teal,
        NamedColor::RebeccaPurple,
        NamedColor::Brown,
        NamedColor::YellowGreen,
        NamedColor::DarkGreen,
        NamedColor::Indigo,
        NamedColor::OliveDrab,
    ];
    let linecol: Vec<NamedColor> = (0..lines_number).map(|l| cols[l]).collect();
    let forecol = Rgb::new(0, 0, 0);
    let gridcol = Rgb::new(220, 220, 220);
    let transp = NamedColor::Transparent;
    let thick: usize = 3;
    let medium: usize = 4;
    let _thin: usize = 2;
    let msize: usize = 10;
    let fsz_title: usize = 35;
    let fsz_legend: usize = 35;
    let fsz_ticks: usize = 30;
    let fsz_axes: usize = 35;

    let mut traces = Vec::new();

    for l in 0..lines_number {
        match plot_par.line_or_points[l] {
            LineOrPoints::Line => {
                traces.push(
                    Scatter::new(x[l].clone(), y[l].clone())
                        .name(&plot_par.legends[l])
                        .mode(Mode::Lines)
                        .line(Line::new().color(linecol[l]).width(medium as f64)),
                )
            },
            LineOrPoints::Points => {
                traces.push(
                    Scatter::new(x[l].clone(), y[l].clone())
                        .name(&plot_par.legends[l])
                        .mode(Mode::Markers)
                        .marker(Marker::new().size(msize).color(linecol[l]).symbol(MarkerSymbol::Circle)),
                )
            },
            LineOrPoints::LineAndPoints => {
                traces.push(
                    Scatter::new(x[l].clone(), y[l].clone())
                        .name(&plot_par.legends[l])
                        .mode(Mode::LinesMarkers)
                        .line(Line::new().color(linecol[l]).width(medium as f64))
                        .marker(Marker::new().size(msize).symbol(MarkerSymbol::Circle)),
                )
            },
        }
    }

    let title = Title::new(&plot_par.title)
        .font(Font::new().size(fsz_title).family("Serif").color(forecol));

    let legend_bottom_right = Legend::new()
        .x(0.99)
        .x_anchor(Anchor::Right)
        .y(0.01)
        .y_anchor(Anchor::Bottom);

    let legend_top_right = Legend::new()
        .x(0.99)
        .x_anchor(Anchor::Right)
        .y(0.99)
        .y_anchor(Anchor::Top);

    let legend_bottom_left = Legend::new()
        .x(0.01)
        .x_anchor(Anchor::Left)
        .y(0.01)
        .y_anchor(Anchor::Bottom);

    let legend_top_left = Legend::new()
        .x(0.01)
        .x_anchor(Anchor::Left)
        .y(0.99)
        .y_anchor(Anchor::Top);

    let legend = match plot_par.legend_al {
        LegendAl::BottomLeft => legend_bottom_left,
        LegendAl::BottomRight => legend_bottom_right,
        LegendAl::TopLeft => legend_top_left,
        LegendAl::TopRight => legend_top_right,
    }.font(Font::new().size(fsz_legend).color(forecol).family("Serif"))
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

    for l in 0..lines_number {
        let trace = traces[l].clone();
        plot.add_trace(trace);
    }
    plot.set_layout(layout);

    plot.write_image(&plot_par.flnm, ImageFormat::PDF, 1280, 960, 1.0);
    plot.write_image(&plot_par.flnm, ImageFormat::PNG, 1280, 960, 1.0);
}