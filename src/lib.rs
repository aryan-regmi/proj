use plotly::common::{
    ColorScale, ColorScalePalette, DashType, Fill, Font, Line, LineShape, Marker, Mode, Title,
};
use plotly::layout::{Axis, BarMode, Layout, Legend, TicksDirection};
use plotly::{Bar, NamedColor, Plot, Rgb, Rgba, Scatter};

// TODO: Rewrite doc comments to look better on cargo docs


// TODO: Add wind + subplots to show wind displacement

/// Converts wind vector to NS(north-south) and EW(east-west) components
fn wind(wind_mag: f64, wind_dir: Degrees) -> (Wind, Wind) {
    let mut direction = wind_dir.0;
    if direction > 360.0 {
        direction = direction % 360.0;
    }
    let north_south = -wind_mag * direction.to_radians().sin();
    let east_west = -wind_mag * direction.to_radians().cos(); 
    
    (Wind::NS(north_south), Wind::EW(east_west))
}

// North = +Y axis
// South = -Y axis
// East = +X axis
// West = -X axis
#[derive(Debug, PartialEq)]
enum Wind {
    NS(f64), // North-South direction (Y axis)wnd
    EW(f64), // East-West direction (X axis)
}

struct Degrees(f64);
impl Degrees {
    fn to_rad(&self) -> Radians {
        Radians(self.0.to_radians())
    }
}

struct Radians(f64);
impl Radians {
    fn to_deg(&self) -> Degrees {
        Degrees(self.0.to_degrees())
    }
}

/// 2D Position Vector
#[allow(dead_code)]
pub struct Position {
    pub x_pos: f64,
    pub y_pos: f64,
}

/// 2D Velocity Vector
#[allow(dead_code)]
pub struct Velocity {
    pub x_vel: f64,
    pub y_vel: f64,
}

/// Defines a Symmetric projectile
///
/// Fields:     mass = Mass of ball
///             cd = Coefficient of Drag for ball
///             area = Cross-sectional area where aerodynamic forces are acting on
pub struct Ball {
    mass: f64,
    cd: f64,
    area: f64,
}
impl Ball {
    fn new(mass: f64, cd: f64, area: f64) -> Self {
        Ball {mass, cd, area}
    }
}

/// Defines an Asymmetric projectile
///
/// Fields:     mass = Mass of ball
///             cd = Coefficient of Drag for ball (x and y directions)
///             area = Cross-sectional area (x and y directions)
pub struct AsymProj {
    mass: f64,
    cd: (f64, f64),
    area: (f64, f64),
}
impl AsymProj {
    fn new(mass: f64, cd: (f64,f64), area: (f64,f64)) -> Self {
        AsymProj {mass, cd, area}
    }
}

pub trait Projectile {
    fn get_mass(&self) -> f64;
    fn is_symmetric(&self) -> bool;
    fn get_cd(&self) -> (f64, f64);
    fn get_area(&self) -> (f64, f64);

    fn trajectory<T: Projectile>(
        projectile: &T,
        init_pos: Position,
        init_vel: Velocity,
        rho: f64,
        g: f64,
        num_iter: usize,
        step_size: f64,
    ) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, usize) {
        // Extract Initial Postions
        let x_curr = init_pos.x_pos;
        let y_curr = init_pos.y_pos;
    
        // Extract Initial Velocities
        let vx_curr = init_vel.x_vel;
        let vy_curr = init_vel.y_vel;
    
        // Extract Projectile Properties
        let mass = projectile.get_mass();
        let (cd_x, cd_y) = projectile.get_cd();
        let (area_x, area_y) = projectile.get_area();
    
        // Initialize Return Values (Set initial values)
        let mut x_vec = vec![0.; num_iter + 1];
        x_vec[0] = x_curr;
        let mut y_vec = vec![0.; num_iter + 1];
        y_vec[0] = y_curr;
        let mut vx = vec![0.; num_iter + 1];
        vx[0] = vx_curr;
        let mut vy = vec![0.; num_iter + 1];
        vy[0] = vy_curr;
        let mut idx = 0;
    
        for i in 1..=num_iter {
            // Drag Force Calculations
            let x_drag = -0.5 * rho * vx_curr * vx_curr * area_x * cd_x * step_size;
            let y_drag = -0.5 * rho * vy_curr * vy_curr * area_y * cd_y * step_size;
    
            // Compute Accelerations
            let ax = (x_drag / mass) * step_size;
            let ay = (-g + (y_drag / mass)) * step_size;
    
            // Update Velocities
            vx[i] = vx[i - 1] + ax;
            vy[i] = vy[i - 1] + ay;
    
            // Update Positions
            x_vec[i] = x_vec[i - 1] + vx[i] * step_size;
            y_vec[i] = y_vec[i - 1] + vy[i] * step_size;
    
            if y_vec[i] < 0.0 {
                idx = i;
                break;
            }
        }
    
        (x_vec, y_vec, vx, vy, idx)
    }

    fn plot_traj(xvec: Vec<f64>, yvec: Vec<f64>, opts: Vec<PlotOpts>) {
        // Default values for plot options
        let (mut xmin, mut xmax, mut ymin, mut ymax) = (0.0, 10.0, 0.0, 10.0);
        let (mut xlabel, mut ylabel) = ("X".to_owned(), "Y".to_owned());
        let mut line_color = NamedColor::Red;
        let mut line_size = 1.5;

        // Read in plot options
        for opt in opts {
            match opt {
                PlotOpts::XMin(val) => xmin = val,
                PlotOpts::XMax(val) => xmax = val,
                PlotOpts::YMin(val) => ymin = val,
                PlotOpts::YMax(val) => ymax = val,
                PlotOpts::XLabel(val) => xlabel = val,
                PlotOpts::YLabel(val) => ylabel = val,
                PlotOpts::LineColor(val) => line_color = val,
                PlotOpts::LineSize(val) => line_size = val,
            }
        }

        let trace = Scatter::new(xvec, yvec)
            .mode(Mode::Lines)
            .line(Line::new().color(line_color).width(line_size));
        let mut plot = Plot::new();
        let layout = Layout::new()
            .title(Title::new("Projectile Trajectory"))
            .x_axis(
                Axis::new()
                    .title(Title::new(&xlabel[..]))
                    .range(vec![xmin, xmax]),
            )
            .y_axis(
                Axis::new()
                    .title(Title::new(&ylabel[..]))
                    .range(vec![ymin, ymax]),
            );

        plot.set_layout(layout);
        plot.add_trace(trace);
        plot.show();
    } 

}

impl Projectile for Ball {
    fn get_mass(&self) -> f64 {
        self.mass
    }
    fn is_symmetric(&self) -> bool {
        true
    }
    fn get_cd(&self) -> (f64, f64) {
        (self.cd, self.cd)
    }
    fn get_area(&self) -> (f64, f64) {
        (self.area, self.area)
    }
}

impl Projectile for AsymProj {
    fn get_mass(&self) -> f64 {
        self.mass
    }
    fn is_symmetric(&self) -> bool {
        false
    }
    fn get_cd(&self) -> (f64, f64) {
        self.cd
    }
    fn get_area(&self) -> (f64, f64) {
        self.area
    }
}

pub enum PlotOpts {
    XMin(f64),
    XMax(f64),
    YMin(f64),
    YMax(f64),
    XLabel(String),
    YLabel(String),
    LineColor(NamedColor),
    LineSize(f64),
}


// TODO: Write tests for AsymProj
// TODO: Integration tests (what it looks like from outside)?
#[cfg(test)]
mod test {
    use super::PlotOpts::*;
    use super::*;
    use speculate::speculate;

    speculate! {

        describe "Projectile Tests" {

            const MASS: f64 = 0.5;
            const CD: f64 = 0.0;
            const CD2: f64 = 0.5;
            const AREA: f64 = 1.0;
            const RHO: f64 = 1.225;
            const G: f64 = 9.8;
            const N: usize = 3000;
            const H: f64 = 0.01;
            const POS: Position = Position {x_pos: 0.0, y_pos: 0.0};
            const VEL: Velocity = Velocity {x_vel: 10.0, y_vel: 10.0};

            /// Rounds given value to specified number of digits
            ///
            /// Parameters:     value = Floating point number to round
            ///                 num_digits = Number of digits to round to (Must be whole number)
            ///
            /// Returns:        Given value with specified number of digits
            ///
            /// *Should only be used to round down the number of digits*
            fn round_dec(value: f64,num_digits: f64) -> f64 {
                (value * 10_f64.powf(num_digits)).round()/(10_f64.powf(num_digits))
            }

            /// Finds max value of a Vec<f64>
            ///
            /// Parameters:     v = A Vector of f64 values to find the max of
            ///
            /// Returns:        The maximum value in the given vector
            fn maxVec(v: Vec<f64>) -> f64 {
                v.iter().cloned().fold(0./0., f64::max)
            }

            before {
                let _ball = Ball::new(MASS, CD, AREA);
                let _ball2 = Ball::new(MASS, CD2, AREA);
            }

            it "implements Projectile trait" {
                assert_eq!(_ball.get_mass(), MASS);
                assert_eq!(_ball.get_cd().0, CD);
                assert_eq!(_ball.get_area().0, AREA);
            }

            it "calculates correct trajectory" {
                let (x,y,_,_,idx) = Ball::trajectory(&_ball, POS, VEL, RHO, G, N, H);
                // Ball::plot_traj(x.clone(), y.clone(), vec![
                //     XMax(20.6),
                //     YMax(6.0),
                //     XLabel("TST".to_owned()),
                //     LineColor(NamedColor::SeaGreen),
                //     LineSize(3.0),
                // ]);
                assert_eq!(round_dec(x[idx],1.), 20.4); // Max Range
                assert_eq!(round_dec(maxVec(y), 1.), 5.1); // Max Height
            }

            it "can factor in air resistance" {
                let (x,y,_,_,idx) = Ball::trajectory(&_ball2, POS, VEL, RHO, G, N, H);
                assert_eq!(round_dec(x[idx],1.), 18.1); // Max Range
                assert_eq!(round_dec(maxVec(y), 1.), 4.8); // Max Height
            }

            it "can convert wind values to components" {
                let wnd = wind(10., Degrees(765.));
                let (mut ns, mut ew) = (wnd.0, wnd.1);
                match ns {
                    Wind::NS(val) => { ns = Wind::NS(round_dec(val,2.)) }
                    Wind::EW(_) => { eprintln!("ERROR: ns should not be Wind::EW") }
                }
                match ew {
                    Wind::EW(val) => { ew = Wind::EW(round_dec(val,2.)) }
                    Wind::NS(_) => { eprintln!("ERROR: ew should not be Wind::NS") }
                }
                assert_eq!((ns, ew), (Wind::NS(-7.07), Wind::EW(-7.07)));
            }
        }

    }
}
