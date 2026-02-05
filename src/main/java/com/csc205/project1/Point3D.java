package com.csc205.project1;

import java.util.logging.Logger;
import java.util.logging.Level;

/**
 * Represents a point in three-dimensional Cartesian space.
 * 
 * This class provides a comprehensive implementation of a 3D point with support for
 * common geometric operations including distance calculations, rotations, translations,
 * and vector operations. The class follows the immutability pattern for thread safety
 * and predictable behavior.
 * 
 * <p>Design Patterns and Principles:</p>
 * <ul>
 *   <li><strong>Immutability Pattern:</strong> All fields are final and the class provides
 *       no setters. Operations return new instances rather than modifying the current object.
 *       This ensures thread safety and prevents unexpected side effects.</li>
 *   <li><strong>Value Object Pattern:</strong> This class represents a mathematical concept
 *       (a point in 3D space) and uses value semantics. Equality is based on the values
 *       of coordinates rather than object identity.</li>
 *   <li><strong>Factory Method Pattern:</strong> Static factory methods like {@code origin()}
 *       provide clear, self-documenting ways to create common point instances.</li>
 *   <li><strong>Single Responsibility Principle:</strong> The class has one clear responsibility:
 *       representing and manipulating 3D points.</li>
 * </ul>
 * 
 * <p>Data Structures & Algorithms Foundations:</p>
 * <ul>
 *   <li><strong>Euclidean Distance Algorithm:</strong> Implements the fundamental distance
 *       formula √((x₂-x₁)² + (y₂-y₁)² + (z₂-z₁)²) with O(1) time complexity.</li>
 *   <li><strong>Rotation Matrices:</strong> Uses linear algebra matrix multiplication for
 *       rotations, demonstrating how mathematical structures map to algorithmic operations.</li>
 *   <li><strong>Vector Operations:</strong> Implements dot product and cross product,
 *       fundamental operations in computational geometry with O(1) complexity.</li>
 *   <li><strong>Normalization:</strong> Demonstrates the concept of unit vectors and
 *       magnitude calculation, essential for many graphics and physics algorithms.</li>
 * </ul>
 * 
 * @author Your Name
 * @version 1.0
 * @since 1.0
 */
public class Point3D {
    
    private static final Logger logger = Logger.getLogger(Point3D.class.getName());
    
    private final double x;
    private final double y;
    private final double z;
    
    // Tolerance for floating-point comparisons
    private static final double EPSILON = 1e-10;
    
    /**
     * Constructs a new Point3D with the specified coordinates.
     * 
     * This is the primary constructor for creating points in 3D space. It uses
     * double precision floating-point numbers to represent coordinates, providing
     * sufficient precision for most scientific and engineering applications.
     * 
     * <p>The constructor validates that none of the coordinates are NaN (Not a Number)
     * or infinite, which would indicate invalid mathematical operations upstream.</p>
     * 
     * <p><strong>Design Note:</strong> This constructor is public to allow flexible
     * point creation, though factory methods are provided for common cases.</p>
     * 
     * @param x the x-coordinate (horizontal axis)
     * @param y the y-coordinate (vertical axis in many coordinate systems)
     * @param z the z-coordinate (depth axis)
     * @throws IllegalArgumentException if any coordinate is NaN or infinite
     */
    public Point3D(double x, double y, double z) {
        if (Double.isNaN(x) || Double.isNaN(y) || Double.isNaN(z)) {
            logger.severe("Attempted to create Point3D with NaN coordinates: " +
                         "x=" + x + ", y=" + y + ", z=" + z);
            throw new IllegalArgumentException("Coordinates cannot be NaN");
        }
        
        if (Double.isInfinite(x) || Double.isInfinite(y) || Double.isInfinite(z)) {
            logger.severe("Attempted to create Point3D with infinite coordinates: " +
                         "x=" + x + ", y=" + y + ", z=" + z);
            throw new IllegalArgumentException("Coordinates cannot be infinite");
        }
        
        this.x = x;
        this.y = y;
        this.z = z;
        
        logger.info("Created new Point3D: " + this);
    }
    
    /**
     * Factory method that creates a Point3D at the origin (0, 0, 0).
     * 
     * This is a convenience factory method that provides a clear, self-documenting
     * way to obtain a point at the coordinate system's origin. Factory methods like
     * this improve code readability and can be optimized (e.g., by caching frequently
     * used instances).
     * 
     * <p><strong>Design Pattern:</strong> Factory Method - provides an alternative
     * to direct constructor calls with more semantic meaning.</p>
     * 
     * @return a new Point3D at coordinates (0, 0, 0)
     */
    public static Point3D origin() {
        logger.info("Creating origin point");
        return new Point3D(0, 0, 0);
    }
    
    /**
     * Calculates the Euclidean distance from this point to another point.
     * 
     * This method implements the fundamental Euclidean distance formula in 3D space:
     * <pre>
     * distance = √((x₂-x₁)² + (y₂-y₁)² + (z₂-z₁)²)
     * </pre>
     * 
     * <p>The algorithm has O(1) time complexity and O(1) space complexity, as it
     * performs a fixed number of arithmetic operations regardless of input size.</p>
     * 
     * <p><strong>Algorithm Foundation:</strong> This is the L2 norm (Euclidean norm)
     * in 3-dimensional space, fundamental to many algorithms including:</p>
     * <ul>
     *   <li>K-nearest neighbors (KNN) classification</li>
     *   <li>K-means clustering</li>
     *   <li>Collision detection in physics engines</li>
     *   <li>Pathfinding algorithms like A*</li>
     * </ul>
     * 
     * @param other the point to calculate distance to; must not be null
     * @return the Euclidean distance between this point and the other point
     * @throws IllegalArgumentException if other is null
     */
    public double distanceTo(Point3D other) {
        if (other == null) {
            logger.severe("Attempted to calculate distance to null point");
            throw new IllegalArgumentException("Cannot calculate distance to null point");
        }
        
        double dx = this.x - other.x;
        double dy = this.y - other.y;
        double dz = this.z - other.z;
        
        double distance = Math.sqrt(dx * dx + dy * dy + dz * dz);
        
        logger.info("Calculated distance from " + this + " to " + other + ": " + distance);
        
        return distance;
    }
    
    /**
     * Calculates the squared distance from this point to another point.
     * 
     * This method computes the square of the Euclidean distance without taking
     * the square root. This is useful for performance-critical applications where
     * only relative distances matter, as it avoids the expensive square root operation.
     * 
     * <p><strong>Performance Optimization:</strong> When comparing distances to find
     * the nearest point, using squared distance is more efficient:</p>
     * <ul>
     *   <li>Square root is computationally expensive (typically 10-20x slower than multiplication)</li>
     *   <li>For comparison purposes: if d₁ &lt; d₂, then d₁² &lt; d₂²</li>
     *   <li>Common in algorithms like KNN, collision detection, and spatial partitioning</li>
     * </ul>
     * 
     * <p><strong>Complexity:</strong> O(1) time, O(1) space</p>
     * 
     * @param other the point to calculate squared distance to; must not be null
     * @return the squared Euclidean distance between this point and the other point
     * @throws IllegalArgumentException if other is null
     */
    public double distanceSquaredTo(Point3D other) {
        if (other == null) {
            logger.severe("Attempted to calculate squared distance to null point");
            throw new IllegalArgumentException("Cannot calculate distance to null point");
        }
        
        double dx = this.x - other.x;
        double dy = this.y - other.y;
        double dz = this.z - other.z;
        
        double distanceSquared = dx * dx + dy * dy + dz * dz;
        
        logger.info("Calculated squared distance from " + this + " to " + other + ": " + distanceSquared);
        
        return distanceSquared;
    }
    
    /**
     * Rotates this point around the X-axis by the specified angle.
     * 
     * This method performs a rotation transformation using the rotation matrix for
     * the X-axis. The rotation follows the right-hand rule: positive angles rotate
     * counterclockwise when looking along the positive X-axis toward the origin.
     * 
     * <p>The rotation matrix for X-axis rotation is:</p>
     * <pre>
     * | 1    0         0      |
     * | 0  cos(θ)  -sin(θ)    |
     * | 0  sin(θ)   cos(θ)    |
     * </pre>
     * 
     * <p><strong>Algorithm Foundation:</strong> This implements matrix-vector
     * multiplication, a fundamental operation in linear algebra with applications in:</p>
     * <ul>
     *   <li>3D graphics and game engines</li>
     *   <li>Robotics and kinematics</li>
     *   <li>Computer vision and image processing</li>
     *   <li>Machine learning transformations</li>
     * </ul>
     * 
     * <p><strong>Immutability:</strong> This method returns a new Point3D rather than
     * modifying the current instance, following the immutability pattern.</p>
     * 
     * @param angleRadians the angle to rotate in radians (not degrees)
     * @return a new Point3D representing the rotated point
     */
    public Point3D rotateX(double angleRadians) {
        if (Double.isNaN(angleRadians) || Double.isInfinite(angleRadians)) {
            logger.warning("Invalid rotation angle: " + angleRadians + ". Returning original point.");
            return this;
        }
        
        double cos = Math.cos(angleRadians);
        double sin = Math.sin(angleRadians);
        
        double newY = y * cos - z * sin;
        double newZ = y * sin + z * cos;
        
        logger.info("Rotated point " + this + " around X-axis by " + 
                   Math.toDegrees(angleRadians) + " degrees");
        
        return new Point3D(x, newY, newZ);
    }
    
    /**
     * Rotates this point around the Y-axis by the specified angle.
     * 
     * This method performs a rotation transformation using the rotation matrix for
     * the Y-axis. The rotation follows the right-hand rule: positive angles rotate
     * counterclockwise when looking along the positive Y-axis toward the origin.
     * 
     * <p>The rotation matrix for Y-axis rotation is:</p>
     * <pre>
     * |  cos(θ)  0  sin(θ)  |
     * |    0     1    0     |
     * | -sin(θ)  0  cos(θ)  |
     * </pre>
     * 
     * <p><strong>Note:</strong> The Y-axis rotation matrix has a different structure
     * than X and Z rotations due to the cyclic nature of the coordinate axes.</p>
     * 
     * @param angleRadians the angle to rotate in radians (not degrees)
     * @return a new Point3D representing the rotated point
     */
    public Point3D rotateY(double angleRadians) {
        if (Double.isNaN(angleRadians) || Double.isInfinite(angleRadians)) {
            logger.warning("Invalid rotation angle: " + angleRadians + ". Returning original point.");
            return this;
        }
        
        double cos = Math.cos(angleRadians);
        double sin = Math.sin(angleRadians);
        
        double newX = x * cos + z * sin;
        double newZ = -x * sin + z * cos;
        
        logger.info("Rotated point " + this + " around Y-axis by " + 
                   Math.toDegrees(angleRadians) + " degrees");
        
        return new Point3D(newX, y, newZ);
    }
    
    /**
     * Rotates this point around the Z-axis by the specified angle.
     * 
     * This method performs a rotation transformation using the rotation matrix for
     * the Z-axis. The rotation follows the right-hand rule: positive angles rotate
     * counterclockwise when looking along the positive Z-axis toward the origin.
     * 
     * <p>The rotation matrix for Z-axis rotation is:</p>
     * <pre>
     * | cos(θ)  -sin(θ)  0 |
     * | sin(θ)   cos(θ)  0 |
     * |   0        0     1 |
     * </pre>
     * 
     * <p><strong>Common Use:</strong> Z-axis rotation is frequently used in 2D graphics
     * (treating Z as "up" out of the plane) and in applications where objects rotate
     * around a vertical axis.</p>
     * 
     * @param angleRadians the angle to rotate in radians (not degrees)
     * @return a new Point3D representing the rotated point
     */
    public Point3D rotateZ(double angleRadians) {
        if (Double.isNaN(angleRadians) || Double.isInfinite(angleRadians)) {
            logger.warning("Invalid rotation angle: " + angleRadians + ". Returning original point.");
            return this;
        }
        
        double cos = Math.cos(angleRadians);
        double sin = Math.sin(angleRadians);
        
        double newX = x * cos - y * sin;
        double newY = x * sin + y * cos;
        
        logger.info("Rotated point " + this + " around Z-axis by " + 
                   Math.toDegrees(angleRadians) + " degrees");
        
        return new Point3D(newX, newY, z);
    }
    
    /**
     * Translates this point by the specified offsets.
     * 
     * Translation is a fundamental geometric transformation that moves a point
     * by adding offset values to each coordinate. This is equivalent to vector
     * addition where the point is treated as a position vector from the origin.
     * 
     * <p><strong>Algorithm:</strong> Simple element-wise addition with O(1) complexity.</p>
     * 
     * <p><strong>Applications:</strong></p>
     * <ul>
     *   <li>Moving objects in 3D space</li>
     *   <li>Camera transformations</li>
     *   <li>Coordinate system changes</li>
     *   <li>Physics simulations (applying velocity/displacement)</li>
     * </ul>
     * 
     * @param dx the offset to add to the x-coordinate
     * @param dy the offset to add to the y-coordinate
     * @param dz the offset to add to the z-coordinate
     * @return a new Point3D at the translated position
     */
    public Point3D translate(double dx, double dy, double dz) {
        if (Double.isNaN(dx) || Double.isNaN(dy) || Double.isNaN(dz)) {
            logger.warning("Invalid translation offsets (NaN detected). Returning original point.");
            return this;
        }
        
        logger.info("Translating point " + this + " by offset (" + dx + ", " + dy + ", " + dz + ")");
        
        return new Point3D(x + dx, y + dy, z + dz);
    }
    
    /**
     * Calculates the magnitude (length) of the vector from the origin to this point.
     * 
     * The magnitude represents the distance from the origin (0, 0, 0) to this point.
     * This is also called the L2 norm or Euclidean norm of the position vector.
     * 
     * <p>Formula: magnitude = √(x² + y² + z²)</p>
     * 
     * <p><strong>Data Structure Concept:</strong> This demonstrates treating a point
     * as a vector, showing the duality between geometric positions and algebraic vectors.</p>
     * 
     * <p><strong>Applications:</strong></p>
     * <ul>
     *   <li>Normalizing vectors to unit length</li>
     *   <li>Calculating speed from velocity vectors</li>
     *   <li>Measuring the "size" of vectors in physics simulations</li>
     * </ul>
     * 
     * @return the magnitude of the vector from origin to this point
     */
    public double magnitude() {
        double mag = Math.sqrt(x * x + y * y + z * z);
        logger.info("Calculated magnitude of " + this + ": " + mag);
        return mag;
    }
    
    /**
     * Returns a normalized version of this point (as a unit vector).
     * 
     * Normalization scales a vector to have a magnitude of 1 while preserving its
     * direction. This is accomplished by dividing each component by the vector's magnitude.
     * 
     * <p><strong>Algorithm:</strong></p>
     * <pre>
     * normalized = (x/mag, y/mag, z/mag) where mag = √(x² + y² + z²)
     * </pre>
     * 
     * <p><strong>Special Cases:</strong></p>
     * <ul>
     *   <li>If the point is at the origin (magnitude = 0), normalization is undefined
     *       mathematically. This method returns the original point and logs a warning.</li>
     *   <li>For very small magnitudes near zero, numerical instability may occur.</li>
     * </ul>
     * 
     * <p><strong>Applications:</strong></p>
     * <ul>
     *   <li>Direction vectors in physics and graphics</li>
     *   <li>Surface normals for lighting calculations</li>
     *   <li>Machine learning feature scaling</li>
     *   <li>Coordinate system basis vectors</li>
     * </ul>
     * 
     * @return a new Point3D with magnitude 1 in the same direction, or this point if magnitude is zero
     */
    public Point3D normalize() {
        double mag = magnitude();
        
        if (mag < EPSILON) {
            logger.warning("Cannot normalize point at or near origin: " + this + 
                          ". Returning original point.");
            return this;
        }
        
        Point3D normalized = new Point3D(x / mag, y / mag, z / mag);
        logger.info("Normalized " + this + " to " + normalized);
        return normalized;
    }
    
    /**
     * Calculates the dot product between this point (as a vector) and another point.
     * 
     * The dot product (also called scalar product or inner product) is a fundamental
     * operation in vector algebra that produces a scalar value. It has numerous
     * geometric and algebraic interpretations.
     * 
     * <p>Formula: a · b = a_x * b_x + a_y * b_y + a_z * b_z</p>
     * 
     * <p><strong>Geometric Interpretation:</strong></p>
     * <pre>
     * a · b = |a| * |b| * cos(θ)
     * where θ is the angle between the vectors
     * </pre>
     * 
     * <p><strong>Applications and Properties:</strong></p>
     * <ul>
     *   <li>Determining if vectors are perpendicular: dot product = 0</li>
     *   <li>Calculating angles between vectors</li>
     *   <li>Projecting one vector onto another</li>
     *   <li>Testing if vectors point in similar directions (positive) or opposite (negative)</li>
     *   <li>Used in lighting calculations (Lambertian reflection)</li>
     *   <li>Machine learning: similarity metrics, neural network computations</li>
     * </ul>
     * 
     * <p><strong>Complexity:</strong> O(1) time, O(1) space</p>
     * 
     * @param other the other point/vector to compute dot product with; must not be null
     * @return the dot product of this vector and the other vector
     * @throws IllegalArgumentException if other is null
     */
    public double dotProduct(Point3D other) {
        if (other == null) {
            logger.severe("Attempted to calculate dot product with null point");
            throw new IllegalArgumentException("Cannot calculate dot product with null point");
        }
        
        double result = this.x * other.x + this.y * other.y + this.z * other.z;
        logger.info("Calculated dot product of " + this + " and " + other + ": " + result);
        return result;
    }
    
    /**
     * Calculates the cross product between this point (as a vector) and another point.
     * 
     * The cross product (also called vector product) produces a new vector that is
     * perpendicular to both input vectors. The magnitude of the result equals the
     * area of the parallelogram formed by the two vectors.
     * 
     * <p>Formula: a × b = (a_y*b_z - a_z*b_y, a_z*b_x - a_x*b_z, a_x*b_y - a_y*b_x)</p>
     * 
     * <p><strong>Properties:</strong></p>
     * <ul>
     *   <li>Result is perpendicular to both input vectors</li>
     *   <li>Anti-commutative: a × b = -(b × a)</li>
     *   <li>Magnitude: |a × b| = |a| * |b| * sin(θ)</li>
     *   <li>Direction follows right-hand rule</li>
     *   <li>Result is zero vector if inputs are parallel</li>
     * </ul>
     * 
     * <p><strong>Applications:</strong></p>
     * <ul>
     *   <li>Calculating surface normals from triangle vertices</li>
     *   <li>Determining if a point is on the left or right of a line (in 2D)</li>
     *   <li>Physics: torque and angular momentum calculations</li>
     *   <li>Graphics: backface culling</li>
     *   <li>Robotics: orientation and rotation calculations</li>
     * </ul>
     * 
     * <p><strong>Complexity:</strong> O(1) time, O(1) space</p>
     * 
     * @param other the other point/vector to compute cross product with; must not be null
     * @return a new Point3D representing the cross product vector
     * @throws IllegalArgumentException if other is null
     */
    public Point3D crossProduct(Point3D other) {
        if (other == null) {
            logger.severe("Attempted to calculate cross product with null point");
            throw new IllegalArgumentException("Cannot calculate cross product with null point");
        }
        
        double newX = this.y * other.z - this.z * other.y;
        double newY = this.z * other.x - this.x * other.z;
        double newZ = this.x * other.y - this.y * other.x;
        
        Point3D result = new Point3D(newX, newY, newZ);
        logger.info("Calculated cross product of " + this + " and " + other + ": " + result);
        return result;
    }
    
    /**
     * Scales this point by a scalar factor.
     * 
     * Scaling multiplies each coordinate by a constant factor, effectively moving
     * the point closer to or farther from the origin along the same direction.
     * 
     * <p><strong>Geometric Effect:</strong></p>
     * <ul>
     *   <li>factor &gt; 1: point moves away from origin</li>
     *   <li>0 &lt; factor &lt; 1: point moves toward origin</li>
     *   <li>factor = 0: point moves to origin</li>
     *   <li>factor &lt; 0: point reflects through origin and scales</li>
     * </ul>
     * 
     * @param factor the scaling factor to apply
     * @return a new Point3D scaled by the given factor
     */
    public Point3D scale(double factor) {
        if (Double.isNaN(factor) || Double.isInfinite(factor)) {
            logger.warning("Invalid scaling factor: " + factor + ". Returning original point.");
            return this;
        }
        
        logger.info("Scaling point " + this + " by factor " + factor);
        return new Point3D(x * factor, y * factor, z * factor);
    }
    
    /**
     * Calculates the midpoint between this point and another point.
     * 
     * The midpoint is the average of corresponding coordinates, representing the
     * point exactly halfway between two points in 3D space.
     * 
     * <p>Formula: midpoint = ((x₁+x₂)/2, (y₁+y₂)/2, (z₁+z₂)/2)</p>
     * 
     * <p><strong>Applications:</strong></p>
     * <ul>
     *   <li>Binary space partitioning trees</li>
     *   <li>Subdivision algorithms (mesh smoothing)</li>
     *   <li>Interpolation between positions</li>
     *   <li>Centroid calculations</li>
     * </ul>
     * 
     * @param other the other point; must not be null
     * @return a new Point3D at the midpoint
     * @throws IllegalArgumentException if other is null
     */
    public Point3D midpoint(Point3D other) {
        if (other == null) {
            logger.severe("Attempted to calculate midpoint with null point");
            throw new IllegalArgumentException("Cannot calculate midpoint with null point");
        }
        
        Point3D mid = new Point3D(
            (this.x + other.x) / 2.0,
            (this.y + other.y) / 2.0,
            (this.z + other.z) / 2.0
        );
        
        logger.info("Calculated midpoint between " + this + " and " + other + ": " + mid);
        return mid;
    }
    
    // Getters
    
    /**
     * Returns the x-coordinate of this point.
     * 
     * @return the x-coordinate
     */
    public double getX() {
        return x;
    }
    
    /**
     * Returns the y-coordinate of this point.
     * 
     * @return the y-coordinate
     */
    public double getY() {
        return y;
    }
    
    /**
     * Returns the z-coordinate of this point.
     * 
     * @return the z-coordinate
     */
    public double getZ() {
        return z;
    }
    
    /**
     * Compares this point to another object for equality.
     * 
     * Two points are considered equal if their corresponding coordinates are equal
     * within a small tolerance (EPSILON) to account for floating-point arithmetic
     * imprecision.
     * 
     * <p><strong>Design Pattern:</strong> Value Object equality - based on field
     * values rather than object identity. This is essential for points to be used
     * correctly in collections like HashSet and as HashMap keys.</p>
     * 
     * <p><strong>Floating-Point Comparison:</strong> Direct equality (==) is problematic
     * for doubles due to rounding errors. This implementation uses an epsilon-based
     * comparison to handle this correctly.</p>
     * 
     * @param obj the object to compare with
     * @return true if the objects represent the same point within tolerance
     */
    @Override
    public boolean equals(Object obj) {
        if (this == obj) return true;
        if (obj == null || getClass() != obj.getClass()) return false;
        
        Point3D other = (Point3D) obj;
        
        return Math.abs(this.x - other.x) < EPSILON &&
               Math.abs(this.y - other.y) < EPSILON &&
               Math.abs(this.z - other.z) < EPSILON;
    }
    
    /**
     * Returns a hash code for this point.
     * 
     * The hash code is computed from the coordinate values to ensure that equal
     * points (according to equals()) have equal hash codes, as required by the
     * Object contract.
     * 
     * <p><strong>Data Structure Foundation:</strong> Hash codes are fundamental to
     * hash-based collections (HashMap, HashSet) which provide O(1) average-case
     * lookup time. Proper hash code implementation is critical for performance.</p>
     * 
     * <p><strong>Implementation:</strong> Uses Double.hashCode() which handles
     * the intricacies of floating-point hashing correctly, then combines the
     * individual hash codes using a prime number multiplier to reduce collisions.</p>
     * 
     * @return a hash code value for this point
     */
    @Override
    public int hashCode() {
        int result = 17;
        result = 31 * result + Double.hashCode(x);
        result = 31 * result + Double.hashCode(y);
        result = 31 * result + Double.hashCode(z);
        return result;
    }
    
    /**
     * Returns a string representation of this point.
     * 
     * The format is "Point3D(x, y, z)" with coordinates displayed to reasonable
     * precision. This is useful for debugging and logging.
     * 
     * @return a string representation of this point
     */
    @Override
    public String toString() {
        return String.format("Point3D(%.3f, %.3f, %.3f)", x, y, z);
    }
}