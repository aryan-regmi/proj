
/// 2D Position Vector
struct Position {
    x_pos: f64,
    y_pos: f64,
}

/// 2D Velocity Vector
struct Velocity {
    x_vel: f64,
    y_vel: f64,
}

struct Ball {
    mass: f64,
    cd: f64,
    area: f64,
}

trait Projectile {
    fn get_mass(&self) -> f64;
    fn get_cd(&self) -> f64;
    fn get_area(&self) -> f64;
    fn new(mass: f64, cd: f64, area: f64) -> Self;

}

impl Projectile for Ball {
    fn get_mass(&self) -> f64 {self.mass}
    fn get_cd(&self) -> f64 {self.cd}
    fn get_area(&self) -> f64 {self.area}
    fn new(mass: f64, cd: f64, area: f64) -> Ball {
        Ball {mass, cd, area}
    }
}

fn trajectory<T: Projectile>(
    projectile: &T, 
    init_pos: Position, 
    init_vel: Velocity, 
    rho: f64, 
    g: f64, 
    num_iter: usize, 
    step_size: f64
) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {

    // Extract Initial Postions
    let x_curr = init_pos.x_pos;
    let y_curr = init_pos.y_pos;

    // Extract Initial Velocities
    let vx_curr = init_vel.x_vel;
    let vy_curr = init_vel.y_vel;

    // Extract Projectile Properties
    let mass = projectile.get_mass();
    let drag_coeff = projectile.get_cd();
    let area = projectile.get_area(); 
    
    // Initialize Return Values (Set initial values)
    let mut x_vec = vec![0.; num_iter];
    x_vec[0] = x_curr;
    let mut y_vec = vec![0.; num_iter];
    y_vec[0] = y_curr;
    let mut vx = vec![0.; num_iter];
    vx[0] = vx_curr;
    let mut vy = vec![0.; num_iter];
    vy[0] = vy_curr;

    for i in 1..=num_iter {
        // Drag Force Calculations
        let x_drag = -0.5 * rho * vx_curr * vx_curr * area * drag_coeff * step_size;
        let y_drag = -0.5 * rho * vy_curr * vy_curr * area * drag_coeff * step_size;

        // Compute Accelerations
        let ax = x_drag/mass * step_size;
        let ay = (-g + (y_drag/mass)) * step_size;

        // Update Velocities
        vx[i] = (vx[i-1] + ax) * step_size;
        vy[i] = (vy[i-1] + ay) * step_size;

        if !(y_vec[i] < 0.0) {
            // Update Positions
            x_vec[i] = (x_vec[i-1] + vx[i]) * step_size;
            y_vec[i] = (y_vec[i-1] + vy[i]) * step_size;
        }
        else {
            break; 
        }
    }

    (x_vec, y_vec, vx, vy)
}

#[cfg(test)]
mod test {
    use speculate::speculate;
    use super::*;


    speculate! {
        
        describe "Projectile Tests" {

            const MASS: f64 = 0.5;
            const CD: f64 = 0.0;
            const AREA: f64 = 1.0;

            before {
                let ball = Ball::new(MASS, CD, AREA);
            }

            it "implements Projectile trait" {
                assert_eq!(ball.get_mass(), MASS);
                assert_eq!(ball.get_cd(), CD);
                assert_eq!(ball.get_area(), AREA);
            }

        }

    }
}
