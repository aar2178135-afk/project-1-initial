package com.csc205.project1;

import java.util.logging.Logger;
import java.util.logging.Level;

/**
 * Represents a line segment in three-dimensional Cartesian space.
 * 
 * This class provides a comprehensive implementation of a 3D line segment with support
 * for common geometric operations including length calculations, distance computations,
 * projections, and line-line relationships. A line is defined by two distinct endpoints
 * in 3D space.
 * 
 * <p><strong>Mathematical Representation:</strong></p>
 * <p>A line segment can be represented parametrically as:</p>
 * <pre>
 * L(t) = P₀ + t(P₁ - P₀), where 0 ≤ t ≤ 1
 * </pre>
 * <p>where P₀ is the start point and P₁ is the end point.</p>
 * 
 * <p><strong>Design Patterns and Principles:</strong></p>
 * <ul>
 *   <li><strong>Immutability Pattern:</strong> All fields are final and the class provides
 *       no setters. Operations return new instances rather than modifying the current object.
 *       This ensures thread safety and prevents unexpected side effects in concurrent
 *       environments.</li>
 *   <li><strong>Value Object Pattern:</strong> This class represents a mathematical concept
 *       (a line segment in 3D space) using value semantics. Equality is based on the values
 *       of the endpoints rather than object identity, making instances safe to use as
 *       map keys or in sets.</li>
 *   <li><strong>Composition Pattern:</strong> The Line3D class is composed of two Point3D
 *       objects rather than using raw coordinates. This promotes code reuse and maintains
 *       a clear separation of concerns.</li>
 *   <li><strong>Factory Method Pattern:</strong> Static factory methods like {@code fromPoints()}
 *       provide clear, self-documenting ways to create line instances and can perform
 *       validation or return cached instances.</li>
 *   <li><strong>Defensive Programming:</strong> Extensive validation of inputs with null
 *       checks and degenerate case handling prevents invalid states and provides clear
 *       error messages through logging.</li>
 *   <li><strong>Single Responsibility Principle:</strong> The class has one clear responsibility:
 *       representing and performing geometric operations on 3D line segments.</li>
 * </ul>
 * 
 * <p><strong>Data Structures & Algorithms Foundations:</strong></p>
 * <ul>
 *   <li><strong>Vector Algebra:</strong> Line operations extensively use vector mathematics,
 *       demonstrating how abstract mathematical structures map to concrete algorithmic
 *       implementations. Direction vectors, dot products, and cross products form the
 *       foundation for distance and projection calculations.</li>
 *   <li><strong>Parametric Representation:</strong> Lines are represented parametrically,
 *       enabling efficient computation of points along the line and nearest point calculations.
 *       This is fundamental to many computational geometry algorithms.</li>
 *   <li><strong>Distance Algorithms:</strong> Implements the skew line distance algorithm,
 *       which uses vector projections to find the shortest distance between two non-parallel
 *       lines in 3D space. This demonstrates the application of linear algebra to geometric
 *       problems with O(1) complexity.</li>
 *   <li><strong>Numerical Stability:</strong> Uses epsilon-based comparisons for handling
 *       floating-point arithmetic issues, a critical consideration in geometric algorithms
 *       to avoid numerical instability and incorrect results.</li>
 *   <li><strong>Geometric Primitives:</strong> Serves as a fundamental building block for
 *       more complex geometric structures like meshes, ray tracers, and collision detection
 *       systems. Many advanced algorithms in graphics and physics decompose complex shapes
 *       into line segments for processing.</li>
 *   <li><strong>Complexity Analysis:</strong> All methods are designed with O(1) time and
 *       space complexity, performing a fixed number of operations regardless of input size,
 *       making them suitable for performance-critical applications.</li>
 * </ul>
 * 
 * <p><strong>Applications:</strong></p>
 * <ul>
 *   <li>Computer Graphics: Ray tracing, wireframe rendering, edge detection</li>
 *   <li>Physics Simulation: Collision detection, trajectory calculation</li>
 *   <li>Robotics: Path planning, sensor ray casting</li>
 *   <li>GIS: Route calculation, proximity analysis</li>
 *   <li>CAD: Geometric modeling, constraint solving</li>
 * </ul>
 * 
 * @author Your Name
 * @version 1.0
 * @since 1.0
 */
public class Line3D {
    
    private static final Logger logger = Logger.getLogger(Line3D.class.getName());
    
    private final Point3D start;
    private final Point3D end;
    
    // Tolerance for floating-point comparisons
    private static final double EPSILON = 1e-10;
    
    /**
     * Constructs a new Line3D with the specified start and end points.
     * 
     * This is the primary constructor for creating line segments in 3D space. The line
     * is defined by two distinct points - the start and end points. The order matters
     * for operations that depend on direction (like direction vector calculation).
     * 
     * <p>The constructor validates that:</p>
     * <ul>
     *   <li>Neither point is null</li>
     *   <li>The two points are not identical (which would create a degenerate line)</li>
     * </ul>
     * 
     * <p><strong>Design Note:</strong> By requiring two distinct points, we prevent the
     * creation of invalid geometric objects that could cause division by zero or other
     * mathematical errors in downstream calculations.</p>
     * 
     * <p><strong>Immutability:</strong> The Point3D objects are stored directly as they
     * are themselves immutable, so no defensive copying is necessary.</p>
     * 
     * @param start the starting point of the line segment; must not be null
     * @param end the ending point of the line segment; must not be null and must be
     *            distinct from start
     * @throws IllegalArgumentException if either point is null or if the points are identical
     */
    public Line3D(Point3D start, Point3D end) {
        if (start == null) {
            logger.severe("Attempted to create Line3D with null start point");
            throw new IllegalArgumentException("Start point cannot be null");
        }
        
        if (end == null) {
            logger.severe("Attempted to create Line3D with null end point");
            throw new IllegalArgumentException("End point cannot be null");
        }
        
        if (start.equals(end)) {
            logger.severe("Attempted to create Line3D with identical start and end points: " + start);
            throw new IllegalArgumentException("Start and end points must be distinct");
        }
        
        this.start = start;
        this.end = end;
        
        logger.info("Created new Line3D from " + start + " to " + end);
    }
    
    /**
     * Factory method that creates a Line3D from two points.
     * 
     * This is a convenience factory method that provides a more semantic alternative to
     * the constructor. It delegates to the primary constructor after documenting the
     * creation intent.
     * 
     * <p><strong>Design Pattern:</strong> Factory Method - provides self-documenting
     * code and allows for future optimization (e.g., caching frequently used lines).</p>
     * 
     * @param start the starting point of the line
     * @param end the ending point of the line
     * @return a new Line3D instance
     * @throws IllegalArgumentException if validation fails
     */
    public static Line3D fromPoints(Point3D start, Point3D end) {
        logger.info("Creating Line3D using factory method");
        return new Line3D(start, end);
    }
    
    /**
     * Calculates the length of this line segment.
     * 
     * This method computes the Euclidean distance between the start and end points,
     * which represents the length of the line segment in 3D space. The calculation
     * uses the distance formula:
     * <pre>
     * length = √((x₂-x₁)² + (y₂-y₁)² + (z₂-z₁)²)
     * </pre>
     * 
     * <p><strong>Algorithm:</strong> Delegates to Point3D's distanceTo method, which
     * implements the Euclidean distance algorithm. This demonstrates code reuse and
     * composition.</p>
     * 
     * <p><strong>Complexity:</strong> O(1) time, O(1) space - performs a constant
     * number of arithmetic operations.</p>
     * 
     * <p><strong>Applications:</strong></p>
     * <ul>
     *   <li>Perimeter calculations in mesh processing</li>
     *   <li>Path length calculations in navigation</li>
     *   <li>Wire length estimation in circuit design</li>
     *   <li>Distance-based filtering and sorting</li>
     * </ul>
     * 
     * @return the length of the line segment as a positive double value
     */
    public double length() {
        double len = start.distanceTo(end);
        logger.info("Calculated length of line " + this + ": " + len);
        return len;
    }
    
    /**
     * Calculates the squared length of this line segment.
     * 
     * This method computes the square of the line's length without taking the square root.
     * This is more efficient than length() when you only need to compare lengths or when
     * the actual length value isn't required.
     * 
     * <p><strong>Performance Optimization:</strong> Avoiding the square root operation
     * provides approximately 10-20x performance improvement. This is crucial in algorithms
     * that compare many line segments, such as:</p>
     * <ul>
     *   <li>Finding the shortest edge in a graph or mesh</li>
     *   <li>Collision detection broad phase filtering</li>
     *   <li>Spatial sorting and partitioning</li>
     * </ul>
     * 
     * <p><strong>Mathematical Property:</strong> Since sqrt is monotonically increasing,
     * if L₁² &lt; L₂², then L₁ &lt; L₂, so squared lengths preserve ordering.</p>
     * 
     * <p><strong>Complexity:</strong> O(1) time, O(1) space</p>
     * 
     * @return the squared length of the line segment
     */
    public double lengthSquared() {
        double lenSq = start.distanceSquaredTo(end);
        logger.info("Calculated squared length of line " + this + ": " + lenSq);
        return lenSq;
    }
    
    /**
     * Returns the direction vector of this line segment.
     * 
     * The direction vector points from the start point to the end point and represents
     * the line's orientation in 3D space. It is computed as:
     * <pre>
     * direction = end - start = (end.x - start.x, end.y - start.y, end.z - start.z)
     * </pre>
     * 
     * <p><strong>Geometric Interpretation:</strong> The direction vector indicates both
     * the direction and magnitude (length) of the line segment. For a unit direction
     * vector, see {@link #unitDirection()}.</p>
     * 
     * <p><strong>Applications:</strong></p>
     * <ul>
     *   <li>Ray casting and ray tracing</li>
     *   <li>Calculating perpendicular vectors</li>
     *   <li>Determining line orientation for lighting calculations</li>
     *   <li>Velocity and acceleration vectors in physics</li>
     * </ul>
     * 
     * <p><strong>Complexity:</strong> O(1) time, O(1) space</p>
     * 
     * @return a Point3D representing the direction vector from start to end
     */
    public Point3D direction() {
        Point3D dir = new Point3D(
            end.getX() - start.getX(),
            end.getY() - start.getY(),
            end.getZ() - start.getZ()
        );
        logger.info("Calculated direction vector for line " + this + ": " + dir);
        return dir;
    }
    
    /**
     * Returns the unit direction vector of this line segment.
     * 
     * The unit direction vector is the normalized direction vector with magnitude 1.
     * It preserves the direction of the line while standardizing its length. This is
     * computed by dividing the direction vector by the line's length:
     * <pre>
     * unitDirection = direction / ||direction||
     * </pre>
     * 
     * <p><strong>Normalization:</strong> Converting a vector to unit length is fundamental
     * in many geometric algorithms where only direction matters, not magnitude. The
     * normalization process ensures numerical stability and consistent scaling.</p>
     * 
     * <p><strong>Applications:</strong></p>
     * <ul>
     *   <li>Surface normal calculations in graphics</li>
     *   <li>Consistent force directions in physics simulations</li>
     *   <li>Direction-based interpolation</li>
     *   <li>Angle calculations using dot products</li>
     * </ul>
     * 
     * <p><strong>Edge Cases:</strong> For very short lines (near zero length), this could
     * result in very large values due to division by a small number. The implementation
     * handles this by checking for near-zero length.</p>
     * 
     * <p><strong>Complexity:</strong> O(1) time, O(1) space</p>
     * 
     * @return a Point3D representing the unit direction vector
     */
    public Point3D unitDirection() {
        double len = length();
        
        if (len < EPSILON) {
            logger.warning("Attempted to get unit direction of near-zero length line: " + this);
            return direction(); // Return un-normalized for degenerate case
        }
        
        Point3D unit = new Point3D(
            (end.getX() - start.getX()) / len,
            (end.getY() - start.getY()) / len,
            (end.getZ() - start.getZ()) / len
        );
        
        logger.info("Calculated unit direction vector for line " + this + ": " + unit);
        return unit;
    }
    
    /**
     * Calculates the midpoint of this line segment.
     * 
     * The midpoint is the point exactly halfway between the start and end points,
     * computed as the average of the two endpoints:
     * <pre>
     * midpoint = (start + end) / 2
     * </pre>
     * 
     * <p><strong>Parametric Representation:</strong> This is equivalent to evaluating
     * the parametric line equation at t = 0.5.</p>
     * 
     * <p><strong>Applications:</strong></p>
     * <ul>
     *   <li>Binary space partitioning (BSP) trees</li>
     *   <li>Subdivision surfaces and mesh refinement</li>
     *   <li>Finding centers of mass for physical objects</li>
     *   <li>Interpolation and animation keyframe generation</li>
     * </ul>
     * 
     * <p><strong>Complexity:</strong> O(1) time, O(1) space</p>
     * 
     * @return a Point3D at the midpoint of the line segment
     */
    public Point3D midpoint() {
        Point3D mid = start.midpoint(end);
        logger.info("Calculated midpoint of line " + this + ": " + mid);
        return mid;
    }
    
    /**
     * Evaluates the line at a specific parameter value using the parametric equation.
     * 
     * This method computes a point along the line using the parametric representation:
     * <pre>
     * L(t) = start + t * (end - start)
     * </pre>
     * 
     * <p><strong>Parameter Interpretation:</strong></p>
     * <ul>
     *   <li>t = 0: returns the start point</li>
     *   <li>t = 0.5: returns the midpoint</li>
     *   <li>t = 1: returns the end point</li>
     *   <li>0 &lt; t &lt; 1: returns a point on the line segment</li>
     *   <li>t &lt; 0 or t &gt; 1: returns a point on the infinite line extension</li>
     * </ul>
     * 
     * <p><strong>Applications:</strong></p>
     * <ul>
     *   <li>Animation and interpolation (tweening)</li>
     *   <li>Bézier curve evaluation (for higher-order curves)</li>
     *   <li>Subdivision algorithms</li>
     *   <li>Ray marching in ray tracers</li>
     *   <li>Generating evenly-spaced sample points</li>
     * </ul>
     * 
     * <p><strong>Complexity:</strong> O(1) time, O(1) space</p>
     * 
     * @param t the parameter value (typically 0 to 1 for points on the segment)
     * @return the point on the line corresponding to parameter t
     */
    public Point3D pointAt(double t) {
        if (Double.isNaN(t) || Double.isInfinite(t)) {
            logger.warning("Invalid parameter t = " + t + " for pointAt. Returning start point.");
            return start;
        }
        
        Point3D dir = direction();
        Point3D point = new Point3D(
            start.getX() + t * dir.getX(),
            start.getY() + t * dir.getY(),
            start.getZ() + t * dir.getZ()
        );
        
        logger.info("Evaluated point at t = " + t + " on line " + this + ": " + point);
        return point;
    }
    
    /**
     * Determines if this line is parallel to another line.
     * 
     * Two lines are parallel if their direction vectors are parallel, which occurs when
     * their cross product is the zero vector (or very close to it, within EPSILON tolerance).
     * Mathematically:
     * <pre>
     * lines are parallel ⟺ direction₁ × direction₂ ≈ 0
     * </pre>
     * 
     * <p><strong>Algorithm:</strong> Uses the cross product test because:</p>
     * <ul>
     *   <li>The cross product of parallel vectors is zero</li>
     *   <li>More numerically stable than comparing individual components</li>
     *   <li>Handles edge cases uniformly</li>
     * </ul>
     * 
     * <p><strong>Geometric Significance:</strong> Parallel lines either never intersect
     * (skew or truly parallel) or are coincident. This test is crucial for determining
     * how to compute line-line distance.</p>
     * 
     * <p><strong>Applications:</strong></p>
     * <ul>
     *   <li>Collision detection preprocessing</li>
     *   <li>Determining if lines will intersect</li>
     *   <li>Architectural CAD constraint solving</li>
     *   <li>Railroad track and road design validation</li>
     * </ul>
     * 
     * <p><strong>Complexity:</strong> O(1) time, O(1) space</p>
     * 
     * @param other the other line to check parallelism with; must not be null
     * @return true if the lines are parallel (within tolerance), false otherwise
     * @throws IllegalArgumentException if other is null
     */
    public boolean isParallel(Line3D other) {
        if (other == null) {
            logger.severe("Attempted to check parallelism with null line");
            throw new IllegalArgumentException("Cannot check parallelism with null line");
        }
        
        Point3D dir1 = this.direction();
        Point3D dir2 = other.direction();
        Point3D cross = dir1.crossProduct(dir2);
        
        // Lines are parallel if cross product is zero (or near-zero)
        double crossMagnitude = cross.magnitude();
        boolean parallel = crossMagnitude < EPSILON;
        
        logger.info("Checked parallelism between " + this + " and " + other + ": " + parallel);
        return parallel;
    }
    
    /**
     * Calculates the shortest distance between this line segment and another line segment.
     * 
     * This method implements a comprehensive algorithm for computing the minimum distance
     * between two line segments in 3D space. The algorithm handles all cases:
     * <ul>
     *   <li>Parallel lines (non-intersecting)</li>
     *   <li>Skew lines (non-parallel, non-intersecting)</li>
     *   <li>Intersecting lines (distance = 0)</li>
     * </ul>
     * 
     * <p><strong>Algorithm Overview:</strong></p>
     * <ol>
     *   <li><strong>Parallel Case:</strong> If lines are parallel, the distance is computed
     *       by finding the perpendicular distance from one line to any point on the other.
     *       This uses vector projection.</li>
     *   <li><strong>Skew Lines Case:</strong> For non-parallel lines, the algorithm computes
     *       the parameters (t, s) for the closest points on each line segment using the
     *       formula derived from minimizing the distance function:
     *       <pre>
     *       D(t,s) = ||L₁(t) - L₂(s)||²
     *       </pre>
     *       The solution involves solving a 2x2 linear system obtained from setting the
     *       partial derivatives to zero.</li>
     *   <li><strong>Clamping:</strong> Since we're working with segments (not infinite lines),
     *       the parameters are clamped to [0,1] to ensure the closest points lie on the
     *       actual segments.</li>
     * </ol>
     * 
     * <p><strong>Mathematical Foundation:</strong></p>
     * <p>For two lines L₁(t) = P₁ + t·D₁ and L₂(s) = P₂ + s·D₂, the closest points occur at:</p>
     * <pre>
     * t = ((P₂-P₁)·D₁(D₂·D₂) - (P₂-P₁)·D₂(D₁·D₂)) / (D₁·D₁)(D₂·D₂) - (D₁·D₂)²
     * s = ((P₂-P₁)·D₂(D₁·D₁) - (P₂-P₁)·D₁(D₁·D₂)) / (D₁·D₁)(D₂·D₂) - (D₁·D₂)²
     * </pre>
     * 
     * <p><strong>Numerical Stability:</strong> The implementation handles several edge cases:</p>
     * <ul>
     *   <li>Near-parallel lines (checked with denominator close to zero)</li>
     *   <li>Very short line segments</li>
     *   <li>Endpoints being closest points</li>
     * </ul>
     * 
     * <p><strong>Applications:</strong></p>
     * <ul>
     *   <li>Collision detection in physics engines and games</li>
     *   <li>Proximity queries in CAD systems</li>
     *   <li>Clearance calculations in robotics path planning</li>
     *   <li>Wire routing in PCB design</li>
     *   <li>Molecular modeling (atomic bond distance analysis)</li>
     *   <li>Network topology analysis (minimum spanning trees)</li>
     * </ul>
     * 
     * <p><strong>Complexity:</strong> O(1) time, O(1) space - performs a fixed number
     * of vector operations and arithmetic calculations.</p>
     * 
     * <p><strong>Related Algorithms:</strong></p>
     * <ul>
     *   <li>Gilbert-Johnson-Keerthi (GJK) algorithm for convex shapes</li>
     *   <li>Separating Axis Theorem (SAT) for collision detection</li>
     *   <li>Closest point on triangle algorithms</li>
     * </ul>
     * 
     * @param other the other line segment; must not be null
     * @return the shortest distance between the two line segments
     * @throws IllegalArgumentException if other is null
     */
    public double shortestDistanceTo(Line3D other) {
        if (other == null) {
            logger.severe("Attempted to calculate distance to null line");
            throw new IllegalArgumentException("Cannot calculate distance to null line");
        }
        
        logger.info("Calculating shortest distance between " + this + " and " + other);
        
        // Get direction vectors
        Point3D dir1 = this.direction();
        Point3D dir2 = other.direction();
        
        // Vector from this.start to other.start
        Point3D w = new Point3D(
            this.start.getX() - other.start.getX(),
            this.start.getY() - other.start.getY(),
            this.start.getZ() - other.start.getZ()
        );
        
        // Compute dot products
        double a = dir1.dotProduct(dir1);  // ||dir1||²
        double b = dir1.dotProduct(dir2);
        double c = dir2.dotProduct(dir2);  // ||dir2||²
        double d = dir1.dotProduct(w);
        double e = dir2.dotProduct(w);
        
        // Denominator for the parameter equations
        double denominator = a * c - b * b;
        
        double t = 0.0; // Parameter for this line
        double s = 0.0; // Parameter for other line
        
        // Check if lines are parallel (denominator near zero)
        if (Math.abs(denominator) < EPSILON) {
            logger.info("Lines are parallel or nearly parallel");
            // Lines are parallel - use simple perpendicular distance
            t = 0.0;
            s = (b > c ? d / b : e / c);
        } else {
            // Lines are not parallel - compute parameters for closest points
            t = (b * e - c * d) / denominator;
            s = (a * e - b * d) / denominator;
        }
        
        // Clamp parameters to [0, 1] for line segments
        t = Math.max(0.0, Math.min(1.0, t));
        s = Math.max(0.0, Math.min(1.0, s));
        
        logger.info("Closest point parameters: t = " + t + ", s = " + s);
        
        // Compute the closest points
        Point3D closestOnThis = this.pointAt(t);
        Point3D closestOnOther = other.pointAt(s);
        
        // Calculate distance between closest points
        double distance = closestOnThis.distanceTo(closestOnOther);
        
        logger.info("Shortest distance: " + distance + 
                   " (closest points: " + closestOnThis + " and " + closestOnOther + ")");
        
        return distance;
    }
    
    /**
     * Calculates the distance from a point to this line segment.
     * 
     * This method computes the shortest (perpendicular) distance from a point to the
     * line segment. The algorithm considers three cases:
     * <ol>
     *   <li>The perpendicular from the point to the infinite line falls on the segment</li>
     *   <li>The closest point is the start endpoint</li>
     *   <li>The closest point is the end endpoint</li>
     * </ol>
     * 
     * <p><strong>Algorithm:</strong></p>
     * <p>The method projects the vector from start to the point onto the line direction.
     * The projection parameter t determines which case applies:</p>
     * <pre>
     * t = (point - start) · direction / ||direction||²
     * 
     * if t < 0: closest point is start
     * if t > 1: closest point is end
     * if 0 ≤ t ≤ 1: closest point is at parameter t on segment
     * </pre>
     * 
     * <p><strong>Geometric Interpretation:</strong> This finds the point on the line segment
     * that is closest to the given point, then computes the Euclidean distance.</p>
     * 
     * <p><strong>Applications:</strong></p>
     * <ul>
     *   <li>Point-in-polygon tests (2D extension)</li>
     *   <li>Snap-to-edge functionality in CAD software</li>
     *   <li>Proximity sensors in robotics</li>
     *   <li>Mouse hover detection in graphics applications</li>
     *   <li>Route deviation measurement in navigation</li>
     * </ul>
     * 
     * <p><strong>Complexity:</strong> O(1) time, O(1) space</p>
     * 
     * @param point the point to measure distance from; must not be null
     * @return the shortest distance from the point to the line segment
     * @throws IllegalArgumentException if point is null
     */
    public double distanceToPoint(Point3D point) {
        if (point == null) {
            logger.severe("Attempted to calculate distance to null point");
            throw new IllegalArgumentException("Cannot calculate distance to null point");
        }
        
        Point3D dir = direction();
        Point3D toPoint = new Point3D(
            point.getX() - start.getX(),
            point.getY() - start.getY(),
            point.getZ() - start.getZ()
        );
        
        // Calculate parameter t for closest point on line
        double lenSquared = lengthSquared();
        double t = toPoint.dotProduct(dir) / lenSquared;
        
        // Clamp t to [0, 1] to stay on segment
        t = Math.max(0.0, Math.min(1.0, t));
        
        Point3D closestPoint = pointAt(t);
        double distance = point.distanceTo(closestPoint);
        
        logger.info("Distance from point " + point + " to line " + this + 
                   ": " + distance + " (closest point: " + closestPoint + ")");
        
        return distance;
    }
    
    /**
     * Projects a point onto this line segment.
     * 
     * This method finds the point on the line segment that is closest to the given point.
     * This is equivalent to dropping a perpendicular from the point to the line and
     * finding where it intersects (or the nearest endpoint if the perpendicular doesn't
     * intersect the segment).
     * 
     * <p><strong>Algorithm:</strong> Uses the same projection technique as distanceToPoint
     * but returns the projected point rather than the distance.</p>
     * 
     * <p><strong>Applications:</strong></p>
     * <ul>
     *   <li>Constraining object movement to a path</li>
     *   <li>Shadow projection in graphics</li>
     *   <li>Force projection in physics simulations</li>
     *   <li>Snapping points to lines in CAD tools</li>
     * </ul>
     * 
     * <p><strong>Complexity:</strong> O(1) time, O(1) space</p>
     * 
     * @param point the point to project; must not be null
     * @return the closest point on the line segment to the given point
     * @throws IllegalArgumentException if point is null
     */
    public Point3D projectPoint(Point3D point) {
        if (point == null) {
            logger.severe("Attempted to project null point");
            throw new IllegalArgumentException("Cannot project null point");
        }
        
        Point3D dir = direction();
        Point3D toPoint = new Point3D(
            point.getX() - start.getX(),
            point.getY() - start.getY(),
            point.getZ() - start.getZ()
        );
        
        double lenSquared = lengthSquared();
        double t = toPoint.dotProduct(dir) / lenSquared;
        
        // Clamp to segment
        t = Math.max(0.0, Math.min(1.0, t));
        
        Point3D projection = pointAt(t);
        logger.info("Projected point " + point + " onto line " + this + ": " + projection);
        
        return projection;
    }
    
    /**
     * Returns a reversed version of this line segment.
     * 
     * Creates a new line with the start and end points swapped, effectively reversing
     * the direction. This is useful when direction matters (e.g., for oriented edges
     * in graphs or normal calculations).
     * 
     * <p><strong>Applications:</strong></p>
     * <ul>
     *   <li>Reversing paths or routes</li>
     *   <li>Flipping surface normals</li>
     *   <li>Bidirectional graph edges</li>
     *   <li>Undo/redo operations in CAD</li>
     * </ul>
     * 
     * @return a new Line3D with start and end swapped
     */
    public Line3D reverse() {
        Line3D reversed = new Line3D(end, start);
        logger.info("Reversed line " + this + " to " + reversed);
        return reversed;
    }
    
    // Getters
    
    /**
     * Returns the starting point of this line segment.
     * 
     * @return the start point
     */
    public Point3D getStart() {
        return start;
    }
    
    /**
     * Returns the ending point of this line segment.
     * 
     * @return the end point
     */
    public Point3D getEnd() {
        return end;
    }
    
    /**
     * Compares this line to another object for equality.
     * 
     * Two lines are considered equal if both their start and end points are equal.
     * The comparison uses Point3D's equals method, which includes epsilon-based
     * tolerance for floating-point comparison.
     * 
     * <p><strong>Design Pattern:</strong> Value Object equality - based on field
     * values rather than object identity. This allows lines to be used correctly
     * in collections.</p>
     * 
     * <p><strong>Note:</strong> This implementation considers Line(A, B) different
     * from Line(B, A) since direction matters. If direction-independent equality
     * is needed, additional logic would be required.</p>
     * 
     * @param obj the object to compare with
     * @return true if the objects represent the same line segment
     */
    @Override
    public boolean equals(Object obj) {
        if (this == obj) return true;
        if (obj == null || getClass() != obj.getClass()) return false;
        
        Line3D other = (Line3D) obj;
        
        boolean equal = this.start.equals(other.start) && this.end.equals(other.end);
        
        if (equal) {
            logger.info("Lines are equal: " + this + " and " + other);
        }
        
        return equal;
    }
    
    /**
     * Returns a hash code for this line.
     * 
     * The hash code is computed from both endpoint values to ensure that equal
     * lines have equal hash codes, as required by the Object contract.
     * 
     * <p><strong>Data Structure Foundation:</strong> Proper hash code implementation
     * is critical for using lines in HashSet and as HashMap keys with O(1) average
     * lookup time.</p>
     * 
     * @return a hash code value for this line
     */
    @Override
    public int hashCode() {
        int result = 17;
        result = 31 * result + start.hashCode();
        result = 31 * result + end.hashCode();
        return result;
    }
    
    /**
     * Returns a string representation of this line.
     * 
     * The format is "Line3D[start -> end]" which clearly shows the line's
     * direction and endpoints. This is useful for debugging and logging.
     * 
     * @return a string representation of this line
     */
    @Override
    public String toString() {
        return String.format("Line3D[%s -> %s]", start, end);
    }
}