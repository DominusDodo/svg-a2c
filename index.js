const TAU = Math.PI * 2
const TO_RADIANS = TAU / 360

/**
 * Calculate the angle between two vectors
 * @returns {Number} The angle between the two vectors
 * @param {Number} ux x-value of first vector
 * @param {Number} uy y-value of first vector
 * @param {Number} vx x-value of second vector
 * @param {Number} vy y-value of second vector
 */
function vectorAngle(ux, uy, vx, vy) {
  const sign = Math.sign(ux * vy - uy * vx)
  let dot  = ux * vx + uy * vy

  // Scale to unit vectors
  dot /= Math.sqrt(ux * ux + uy * uy) * Math.sqrt(vx * vx + vy * vy)

  // Adjust for rounding errors
  if (dot >  1.0) dot =  1.0
  if (dot < -1.0) dot = -1.0

  return sign * Math.acos(dot)
}

/**
 * (F.6.5) Conversion from endpoint to center parameterization
 * 
 * Given the following variables:
 *    x1, y1, x2, y2, fa, fs, rx, ry, phi  
 * the task is to find:
 *    cx, cy, theta1, delta_theta
 * 
 * @see https://www.w3.org/TR/SVG11/implnote.html#ArcConversionEndpointToCenter
 * @returns {[Number, Number, Number, Number]} [cx, cy, theta1, delta_theta]
 * @param {Number} x1 x-value of the start point
 * @param {Number} y1 y-value of the start point
 * @param {Number} x2 x-value of the final point
 * @param {Number} y2 y-value of the final point
 * @param {Number} fa large-arc flag
 * @param {Number} fs sweep flag
 * @param {Number} rx x-radius
 * @param {Number} ry y-radius
 * @param {Number} phi angle to the x-axis
 */
function f_6_5(x1, y1, x2, y2, fa, fs, rx, ry, phi) {
  // Compute trigonometric values
  const sin_phi = Math.sin(phi * TAU / 360)
  const cos_phi = Math.cos(phi * TAU / 360)

  /*
   * (F.6.5.1) Step 1
   *
   * Compute (x1', y1')
   *
   * (x1')   (  cos[phi]  sin[phi] ) . ( [x1-x2] / 2 )
   * (x2') = ( -sin[phi]  cos[phi] )   ( [y1-y2] / 2 )
   */
  const x1_m_x2_d_2 = (x1 - x2) / 2
  const y1_m_y2_d_2 = (y1 - y2) / 2
  const x1p =  cos_phi * x1_m_x2_d_2 + sin_phi * y1_m_y2_d_2
  const y1p = -sin_phi * x1_m_x2_d_2 + cos_phi * y1_m_y2_d_2

  /*
   * (F.6.5.2) Step 2
   *
   * Compute (cx', cy')
   * 
   * (cx')       (  [rx*y1'] / ry )  
   * (cy') = k * ( -[ry*x1'] / rx )
   * 
   * k = sign * sqrt(a/b)
   * a = [rx^2 * ry^2] - b
   * b = [rx^2 * y1'^2] - [ry^2 * x1'^2]
   * 
   * sign = { -1 :   if fa = fs
   *        {  1 :   otherwise
   */
  const rx_sq  =  rx * rx
  const ry_sq  =  ry * ry
  const x1p_sq = x1p * x1p
  const y1p_sq = y1p * y1p

  const b = (rx_sq * y1p_sq) + (ry_sq * x1p_sq)
  const a = (rx_sq * ry_sq) - b
  const sign = fa === fs ? -1 : 1

  // Use Math.max to adjust for rounding error
  const k = sign * Math.sqrt(Math.max(0, a / b))

  const cxp = k *  rx * y1p / ry
  const cyp = k * -ry * x1p / rx

  /*
   * (F.6.5.3) Step 3
   *
   * Compute (cx, cy) from (cx', cy')
   * 
   * (x1')   (  cos[phi]  sin[phi] ) . (cx') + ( [x1+x2] / 2 )
   * (x2') = ( -sin[phi]  cos[phi] )   (cx') + ( [y1+y2] / 2 )
   */
  const cx = cos_phi * cxp - sin_phi * cyp + (x1 + x2) / 2
  const cy = sin_phi * cxp + cos_phi * cyp + (y1 + y2) / 2

  /*
   * (F.6.5.5) Step 4
   * (F.6.5.6)
   * 
   * Compute theta1 and delta_theta
   * 
   *          (1)       ( [x1'- cx'] / rx )
   * theta1 = (0) angle ( [y1'- cy'] / ry )
   * 
   *         ( [x1'- cx'] / rx )       ( [-x1'- cx'] / rx )
   * delta = ( [y1'- cy'] / ry ) angle ( [-y1'- cy'] / ry ) mod TAU
   */
  const v1x = ( x1p - cxp) / rx
  const v1y = ( y1p - cyp) / ry
  const v2x = (-x1p - cxp) / rx
  const v2y = (-y1p - cyp) / ry

  const theta1 = vectorAngle(1, 0, v1x, v1y)
  console.log({x1, y1, x2, y2, rx, ry, cos_phi, sin_phi, phi, x1p, y1p, a, b, k, sign,
    x1_m_x2_d_2, y1_m_y2_d_2})
  let delta_theta = vectorAngle(v1x, v1y, v2x, v2y)

  // mod TAU
  if (!fs && delta_theta > 0) delta_theta -= TAU
  if ( fs && delta_theta < 0) delta_theta += TAU

  // Return results
  return [ cx, cy, theta1, delta_theta ]
}

/**
 * (F.6.6) Correction of out-of-range radii
 * 
 * Corrects out-of-range radii as described in (F.6.2) and formalized in (F.6.6)
 * 
 * @see https://www.w3.org/TR/SVG11/implnote.html#ArcCorrectionOutOfRangeRadii
 * @returns {[Number, Number]} [rx, ry]
 * @param {Number} rx x-radius
 * @param {Number} ry y-radius
 * @param {Number} x1p transformed x-value of the start point
 * @param {Number} y1p transformed y-value of the start point
 */
function f_6_6(rx, ry, x1p, y1p) {
  /*
   * (F.6.6) Step 1
   *
   * Ensure radii are non-zero
   * 
   * If rx = 0 or ry = 0, then treat this as a straight line from
   * (x1, y1) to (x2, y2) and stop
   */
  if (rx == 0 || ry == 0) return

  /*
   * (F.6.6.1) Step 2
   *
   * Ensure radii are positive
   * 
   * rx <- |rx|
   * ry <- |ry|
   */
  rx = Math.abs(rx)
  ry = Math.abs(ry)

  /*
   * (F.6.6.2) Step 3
   * (F.6.6.3)
   *
   * Ensure radii are large enough
   * 
   * lambda = (x1'^2 / rx^2) + (y1'^2 / ry^2)
   * 
   * if [lambda > 1]:
   * 
   *    ( rx )                   ( rx )
   *    ( ry ) <- sqrt[lambda] * ( ry )
   */
  const rx_sq  =  rx * rx
  const ry_sq  =  ry * ry
  const x1p_sq = x1p * x1p
  const y1p_sq = y1p * y1p

  const lambda = x1p_sq / rx_sq + y1p_sq / ry_sq

  if (lambda > 1) {
    const sqrt_lambda = Math.sqrt(lambda)
    rx *= sqrt_lambda
    ry *= sqrt_lambda
  }

  // Return corrected values
  return [rx, ry]
}

/**
 * Approximate one unit arc segment with bézier curves
 * @see http://math.stackexchange.com/questions/873224
 * @returns {Number[]} The bezier vertices approximating this unit arc
 * @param {Number} theta1 
 * @param {Number} delta_theta 
 */
function approximateUnitArc(theta1, delta_theta, full) {
  const alpha = 4/3 * Math.tan(delta_theta / 4)

  const x1 = Math.cos(theta1)
  const y1 = Math.sin(theta1)
  const x2 = Math.cos(theta1 + delta_theta)
  const y2 = Math.sin(theta1 + delta_theta)

  const first = full ? [x1, y1] : []

  return [
    ...first,
    x1 - y1 * alpha, y1 + x1 * alpha,
    x2 + y2 * alpha, y2 - x2 * alpha,
    x2, y2
  ]
}

/**
 * Finds an approximation of a given arc using (cubic) bezier curves
 * @returns {Number[][][]} An array of curves (array of 3 or 4 points (array of 2 numbers))
 * @param {Number} x1 x-value of the start point
 * @param {Number} y1 y-value of the start point
 * @param {Number} rx x-radius
 * @param {Number} ry y-radius
 * @param {Number} phi angle to the x-axis
 * @param {Number} fa large-arc flag
 * @param {Number} fs sweep flag
 * @param {Number} x2 x-value of the final point
 * @param {Number} y2 y-value of the final point
 * @param {Boolean} full Whether to return full bezier descriptions
 */
function a2c(x1, y1, rx, ry, phi, fa, fs, x2, y2, full=false) {
  // Calculate trigonometric values
  const sin_phi = Math.sin(phi * TO_RADIANS)
  const cos_phi = Math.cos(phi * TO_RADIANS)

  // Make sure radii are valid
  const x1_m_x2_d_2 = (x1 - x2) / 2
  const y1_m_y2_d_2 = (y1 - y2) / 2
  const x1p =  cos_phi * x1_m_x2_d_2 + sin_phi * y1_m_y2_d_2
  const y1p = -sin_phi * x1_m_x2_d_2 + cos_phi * y1_m_y2_d_2
  
  // Correct out-of-range radii
  const rxry = f_6_6(rx, ry, x1p, y1p)
  // Stop if invalid radii
  if (!rxry) return
  // Extract correct values
  [rx, ry] = rxry

  // Get center parameters (cx, cy, theta1, delta_theta)
  const [cx, cy, theta1, delta_theta] = f_6_5(x1, y1, x2, y2, fa, fs, rx, ry, sin_phi, cos_phi)

  // Split an arc into multiple segments, so that each
  // segment will be less than 90 degrees (TAU / 4)
  const segments = Math.max(Math.ceil(Math.abs(delta_theta) / (TAU / 4)), 1)
  const segmentLength = delta_theta / segments

  // Each 'segment' is a curve
  return Array(segments).fill().map((_, segment) => {
    // Starting angle of this arc approximation
    const theta = theta1 + segmentLength * segment
    // Approximate bezier from arc at theta, using radial length
    const rawCurve = approximateUnitArc(theta, segmentLength, full)
    const curve = []

    // Transform back to original ellipse
    for (let i = 0; i < rawCurve.length; i += 2) {
      // Scale
      const x = rawCurve[i + 0] * rx
      const y = rawCurve[i + 1] * ry
      // Rotate
      const xp = cos_phi * x - sin_phi * y
      const yp = sin_phi * x + cos_phi * y
      // Translate and push
      curve.push([xp + cx, yp + cy])
    }

    // Return transformed curve
    return curve
  })
}

module.exports = a2c
