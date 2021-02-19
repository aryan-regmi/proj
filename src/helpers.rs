/// Rounds given value to specified number of digits
/// 
/// =============PARAMETERS===========
/// value = Floating point number to round
/// num_digits = Number of digits to round to (Must be whole number)
///
/// ==============RETURNS==============
/// Given value with specified number of digits
///
/// *Should only be used to round down the number of digits*
fn round_dec(value: f64,num_digits: f64) -> f64 {
    (value * 10_f64.powf(num_digits)).round()/(10_f64.powf(num_digits))
}

/// Finds max value of a Vec<f64>
///
/// =============PARAMETERS===========
/// v = A Vector of f64 values to find the max of
///
/// ==============RETURNS==============
/// The maximum value in the given vector
fn maxVec(v: Vec<f64>) -> f64 {
    v.iter().cloned().fold(0./0., f64::max)
}
