package com.csc205.project1;

import java.util.logging.Logger;
import java.util.logging.Level;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Represents a cube in three-dimensional Cartesian space.
 * 
 * This class provides a comprehensive implementation of a 3D cube (rectangular
 * parallelepiped) with support for geometric operations including volume calculations,
 * surface area, diagonal measurements, transformations (rotation, translation, scaling),
 * and spatial queries. The cube is defined by eight vertices positioned in 3D space.
 * 
 * <p><strong>Mathematical Representation:</strong></p>
 * <p>A cube is represented by eight vertices that form the corners of a rectangular
 * solid. The vertices follow standard cube vertex numbering with consistent winding.</p>
 * 
 * <p><strong>Design Patterns and Principles:</strong></p>
 * <ul>
 *   <li><strong>Immutability Pattern:</strong> All fields are final, operations return new
 *       instances ensuring thread safety.</li>
 *   <li><strong>Value Object Pattern:</strong> Equality based on vertex positions.</li>
 *   <li><strong>Composition Pattern:</strong> Aggregates Point3D and uses Line3D.</li>
 *   <li><strong>Builder Pattern:</strong> Factory methods for different construction modes.</li>
 *   <li><strong>Collection Pattern:</strong> Returns unmodifiable collections.</li>
 *   <li><strong>Single Responsibility Principle:</strong> Focused on 3D cubes only.</li>
 * </ul>
 * 
 * <p><strong>Data Structures & Algorithms Foundations:</strong></p>
 * <ul>
 *   <li><strong>Spatial Data Structure:</strong> Represents complex 3D objects as collections
 *       of primitives (vertices/edges) - fundamental to mesh representations.</li>
 *   <li><strong>Topological Connectivity:</strong> Maintains graph structure (12 edges, 8
 *       vertices) showing geometry-graph theory connection.</li>
 *   <li><strong>Bounding Volume (AABB):</strong> Axis-aligned bounding boxes enable O(1)
 *       broad-phase collision detection in physics engines.</li>
 *   <li><strong>Transformation Matrices:</strong> Linear algebra for rotations/scaling.</li>
 *   <li><strong>Geometric Calculations:</strong> All O(1) using direct formulas.</li>
 *   <li><strong>Containment Testing:</strong> O(1) point-in-box tests for spatial queries.</li>
 * </ul>
 * 
 * <p><strong>Applications in 3D Graphics:</strong></p>
 * Bounding Volume Hierarchies, collision detection, frustum culling, voxelization,
 * spatial partitioning (octrees), LOD selection, object picking, 3D modeling primitives.
 * 
 * @author Your Name
 * @version 1.0
 * @since 1.0
 */
public class Cube3D {
    
    private static final Logger logger = Logger.getLogger(Cube3D.class.getName());
    
    // Eight vertices defining the cube  
    private final Point3D v0, v1, v2, v3, v4, v5, v6, v7;
    
    // Tolerance for floating-point comparisons
    private static final double EPSILON = 1e-10;
    
    /**
     * Constructs a new Cube3D with the specified eight vertices.
     * 
     * This is the primary constructor that directly accepts all eight vertices of the cube.
     * The vertices must be ordered: vertices 0-3 form front face, vertices 4-7 form back face.
     * 
     * <p><strong>Validation:</strong> Validates that all vertices are non-null but does not
     * enforce valid cube geometry. For guaranteed valid cubes, use the factory methods.</p>
     * 
     * @param v0 front-bottom-left vertex; must not be null
     * @param v1 front-bottom-right vertex; must not be null
     * @param v2 front-top-right vertex; must not be null
     * @param v3 front-top-left vertex; must not be null
     * @param v4 back-bottom-left vertex; must not be null
     * @param v5 back-bottom-right vertex; must not be null
     * @param v6 back-top-right vertex; must not be null
     * @param v7 back-top-left vertex; must not be null
     * @throws IllegalArgumentException if any vertex is null
     */
    public Cube3D(Point3D v0, Point3D v1, Point3D v2, Point3D v3,
                  Point3D v4, Point3D v5, Point3D v6, Point3D v7) {
        
        if (v0 == null || v1 == null || v2 == null || v3 == null ||
            v4 == null || v5 == null || v6 == null || v7 == null) {
            logger.severe("Attempted to create Cube3D with one or more null vertices");
            throw new IllegalArgumentException("All vertices must be non-null");
        }
        
        this.v0 = v0;
        this.v1 = v1;
        this.v2 = v2;
        this.v3 = v3;
        this.v4 = v4;
        this.v5 = v5;
        this.v6 = v6;
        this.v7 = v7;
        
        logger.info("Created new Cube3D with 8 vertices");
    }
    
    /**
     * Factory method to create a cube from center point and dimensions.
     * 
     * This is the most convenient way to create an axis-aligned cube. The cube is
     * constructed by placing vertices at appropriate offsets from the center point.
     * 
     * <p><strong>Design Pattern:</strong> Factory Method - provides a clear, semantic
     * way to construct cubes that aligns with typical 3D thinking (center + size).</p>
     * 
     * <p><strong>Applications:</strong> Creating bounding boxes, defining regions,
     * procedural generation, voxel representations.</p>
     * 
     * @param center the center point of the cube; must not be null
     * @param width the dimension along the X-axis; must be positive
     * @param height the dimension along the Y-axis; must be positive
     * @param depth the dimension along the Z-axis; must be positive
     * @return a new Cube3D instance centered at the specified point
     * @throws IllegalArgumentException if center is null or any dimension is non-positive
     */
    public static Cube3D fromCenterAndDimensions(Point3D center, double width, 
                                                  double height, double depth) {
        if (center == null) {
            logger.severe("Attempted to create cube with null center point");
            throw new IllegalArgumentException("Center point cannot be null");
        }
        
        if (width <= 0 || height <= 0 || depth <= 0) {
            logger.severe("Attempted to create cube with non-positive dimensions");
            throw new IllegalArgumentException("All dimensions must be positive");
        }
        
        double halfW = width / 2.0;
        double halfH = height / 2.0;
        double halfD = depth / 2.0;
        
        double cx = center.getX();
        double cy = center.getY();
        double cz = center.getZ();
        
        Point3D v0 = new Point3D(cx - halfW, cy - halfH, cz - halfD);
        Point3D v1 = new Point3D(cx + halfW, cy - halfH, cz - halfD);
        Point3D v2 = new Point3D(cx + halfW, cy + halfH, cz - halfD);
        Point3D v3 = new Point3D(cx - halfW, cy + halfH, cz - halfD);
        Point3D v4 = new Point3D(cx - halfW, cy - halfH, cz + halfD);
        Point3D v5 = new Point3D(cx + halfW, cy - halfH, cz + halfD);
        Point3D v6 = new Point3D(cx + halfW, cy + halfH, cz + halfD);
        Point3D v7 = new Point3D(cx - halfW, cy + halfH, cz + halfD);
        
        logger.info("Created cube from center with dimensions " + width + "x" + height + "x" + depth);
        
        return new Cube3D(v0, v1, v2, v3, v4, v5, v6, v7);
    }
    
    /**
     * Factory method to create a cube from two opposite corner points.
     * 
     * Given any two opposite corners of an axis-aligned cube, this method constructs
     * the complete cube by computing the remaining vertices. Useful for bounding boxes.
     * 
     * <p><strong>Algorithm:</strong> Identifies min/max X, Y, Z coordinates from the two
     * corners to determine the bounding box, then creates all eight vertices.</p>
     * 
     * <p><strong>Applications:</strong> Computing bounding boxes from point clouds,
     * region selection, defining search spaces for spatial queries.</p>
     * 
     * @param corner1 one corner of the cube; must not be null
     * @param corner2 the opposite corner; must not be null and not equal to corner1
     * @return a new Cube3D instance
     * @throws IllegalArgumentException if either corner is null or corners are identical
     */
    public static Cube3D fromOppositeCorners(Point3D corner1, Point3D corner2) {
        if (corner1 == null || corner2 == null) {
            logger.severe("Attempted to create cube with null corner(s)");
            throw new IllegalArgumentException("Corners cannot be null");
        }
        
        if (corner1.equals(corner2)) {
            logger.severe("Attempted to create cube with identical corners");
            throw new IllegalArgumentException("Corners must be distinct");
        }
        
        double minX = Math.min(corner1.getX(), corner2.getX());
        double maxX = Math.max(corner1.getX(), corner2.getX());
        double minY = Math.min(corner1.getY(), corner2.getY());
        double maxY = Math.max(corner1.getY(), corner2.getY());
        double minZ = Math.min(corner1.getZ(), corner2.getZ());
        double maxZ = Math.max(corner1.getZ(), corner2.getZ());
        
        Point3D v0 = new Point3D(minX, minY, minZ);
        Point3D v1 = new Point3D(maxX, minY, minZ);
        Point3D v2 = new Point3D(maxX, maxY, minZ);
        Point3D v3 = new Point3D(minX, maxY, minZ);
        Point3D v4 = new Point3D(minX, minY, maxZ);
        Point3D v5 = new Point3D(maxX, minY, maxZ);
        Point3D v6 = new Point3D(maxX, maxY, maxZ);
        Point3D v7 = new Point3D(minX, maxY, maxZ);
        
        logger.info("Created cube from opposite corners");
        
        return new Cube3D(v0, v1, v2, v3, v4, v5, v6, v7);
    }
    
    /**
     * Factory method to create a unit cube at the origin.
     * 
     * Creates a 1x1x1 cube centered at the origin (0, 0, 0). Useful for testing,
     * as a base primitive, or as a reference object for scaling/positioning.
     * 
     * @return a new Cube3D instance representing a unit cube at origin
     */
    public static Cube3D unitCube() {
        logger.info("Creating unit cube at origin");
        return fromCenterAndDimensions(Point3D.origin(), 1.0, 1.0, 1.0);
    }
    
    /**
     * Calculates the volume of this cube.
     * 
     * The volume is computed as width × height × depth. For an axis-aligned
     * rectangular parallelepiped, this is a direct O(1) calculation.
     * 
     * <p><strong>Formula:</strong> V = width × height × depth</p>
     * 
     * <p><strong>Applications:</strong> Mass calculations (volume × density), memory
     * estimation for voxels, capacity planning in 3D bin packing, LOD selection.</p>
     * 
     * <p><strong>Complexity:</strong> O(1) time, O(1) space</p>
     * 
     * @return the volume of the cube as a positive double value
     */
    public double volume() {
        double width = v1.distanceTo(v0);
        double height = v3.distanceTo(v0);
        double depth = v4.distanceTo(v0);
        
        double vol = width * height * depth;
        
        logger.info("Calculated volume of cube: " + vol);
        
        return vol;
    }
    
    /**
     * Calculates the surface area of this cube.
     * 
     * The surface area is the sum of all six faces: SA = 2(wh + wd + hd)
     * where w = width, h = height, d = depth.
     * 
     * <p><strong>Applications:</strong> Material estimation for 3D printing, heat
     * transfer calculations, texture mapping memory, collision surface area.</p>
     * 
     * <p><strong>Complexity:</strong> O(1) time, O(1) space</p>
     * 
     * @return the surface area of the cube as a positive double value
     */
    public double surfaceArea() {
        double width = v1.distanceTo(v0);
        double height = v3.distanceTo(v0);
        double depth = v4.distanceTo(v0);
        
        double area = 2.0 * (width * height + width * depth + height * depth);
        
        logger.info("Calculated surface area of cube: " + area);
        
        return area;
    }
    
    /**
     * Calculates the length of the space diagonal of the cube.
     * 
     * The space diagonal (body diagonal) connects two opposite vertices through the
     * cube's interior. For a rectangular parallelepiped: diagonal = √(w² + h² + d²)
     * 
     * <p><strong>Geometric Interpretation:</strong> This is the longest straight line
     * that fits inside the cube.</p>
     * 
     * <p><strong>Applications:</strong> Sphere-box collision (compare radius to half
     * diagonal), camera distance calculations, bounding sphere radius.</p>
     * 
     * <p><strong>Complexity:</strong> O(1) time, O(1) space</p>
     * 
     * @return the length of the space diagonal
     */
    public double spaceDiagonal() {
        double diagonal = v0.distanceTo(v6);
        
        logger.info("Calculated space diagonal of cube: " + diagonal);
        
        return diagonal;
    }
    
    /**
     * Calculates the total perimeter (edge length) of the cube.
     * 
     * A cube has 12 edges: 4 along each principal axis. The total perimeter is
     * perimeter = 4 × (width + height + depth)
     * 
     * <p><strong>Why This Matters:</strong> In wireframe rendering, perimeter represents
     * the total length of line segments to draw. Useful for:</p>
     * <ul>
     *   <li>Estimating wireframe rendering cost</li>
     *   <li>Material length calculations for physical models</li>
     *   <li>Edge traversal algorithms in mesh processing</li>
     * </ul>
     * 
     * <p><strong>Complexity:</strong> O(1) time, O(1) space</p>
     * 
     * @return the total edge length (perimeter) of the cube
     */
    public double perimeter() {
        double width = v1.distanceTo(v0);
        double height = v3.distanceTo(v0);
        double depth = v4.distanceTo(v0);
        
        double perim = 4.0 * (width + height + depth);
        
        logger.info("Calculated perimeter of cube: " + perim);
        
        return perim;
    }
    
    /**
     * Returns the center point of the cube.
     * 
     * The center is computed as the average of all eight vertices.
     * For an axis-aligned cube, this equals the midpoint between opposite vertices.
     * 
     * <p><strong>Algorithm:</strong> center = (v0 + v1 + ... + v7) / 8</p>
     * 
     * <p><strong>Applications:</strong> Object positioning, transformation pivots,
     * centroid-based spatial partitioning, center of mass (uniform density),
     * distance-based sorting.</p>
     * 
     * <p><strong>Complexity:</strong> O(1) time, O(1) space</p>
     * 
     * @return a new Point3D representing the center of the cube
     */
    public Point3D getCenter() {
        double centerX = (v0.getX() + v1.getX() + v2.getX() + v3.getX() +
                         v4.getX() + v5.getX() + v6.getX() + v7.getX()) / 8.0;
        double centerY = (v0.getY() + v1.getY() + v2.getY() + v3.getY() +
                         v4.getY() + v5.getY() + v6.getY() + v7.getY()) / 8.0;
        double centerZ = (v0.getZ() + v1.getZ() + v2.getZ() + v3.getZ() +
                         v4.getZ() + v5.getZ() + v6.getZ() + v7.getZ()) / 8.0;
        
        Point3D center = new Point3D(centerX, centerY, centerZ);
        
        logger.info("Calculated center of cube: " + center);
        
        return center;
    }
    
    /**
     * Returns the dimensions of the cube as an array [width, height, depth].
     * 
     * The dimensions are extracted by measuring distances between appropriate
     * vertices along each principal axis.
     * 
     * @return a 3-element array containing [width, height, depth]
     */
    public double[] getDimensions() {
        double width = v1.distanceTo(v0);
        double height = v3.distanceTo(v0);
        double depth = v4.distanceTo(v0);
        
        return new double[] { width, height, depth };
    }
    
    /**
     * Translates the cube by the specified offset.
     * 
     * Translation moves the entire cube by adding the same offset to all vertices.
     * This is a rigid body transformation preserving size, shape, and orientation.
     * 
     * <p><strong>Algorithm:</strong> Applies translation vector to each of the eight
     * vertices independently, creating a new cube at the translated position.</p>
     * 
     * <p><strong>Transformation Properties:</strong></p>
     * <ul>
     *   <li>Preserves: volume, surface area, edge lengths, angles</li>
     *   <li>Changes: position of all vertices and center point</li>
     *   <li>Composition: T₁ followed by T₂ equals T₁₊₂ (commutative)</li>
     * </ul>
     * 
     * <p><strong>Applications:</strong> Object movement in 3D scenes, camera
     * transformations, physics integration (velocity/displacement), procedural placement.</p>
     * 
     * <p><strong>Complexity:</strong> O(1) time (8 vertex translations), O(1) space</p>
     * 
     * @param dx the offset along the X-axis
     * @param dy the offset along the Y-axis
     * @param dz the offset along the Z-axis
     * @return a new Cube3D at the translated position
     */
    public Cube3D translate(double dx, double dy, double dz) {
        if (Double.isNaN(dx) || Double.isNaN(dy) || Double.isNaN(dz)) {
            logger.warning("Invalid translation offsets (NaN). Returning original cube.");
            return this;
        }
        
        Point3D newV0 = v0.translate(dx, dy, dz);
        Point3D newV1 = v1.translate(dx, dy, dz);
        Point3D newV2 = v2.translate(dx, dy, dz);
        Point3D newV3 = v3.translate(dx, dy, dz);
        Point3D newV4 = v4.translate(dx, dy, dz);
        Point3D newV5 = v5.translate(dx, dy, dz);
        Point3D newV6 = v6.translate(dx, dy, dz);
        Point3D newV7 = v7.translate(dx, dy, dz);
        
        logger.info("Translated cube by offset (" + dx + ", " + dy + ", " + dz + ")");
        
        return new Cube3D(newV0, newV1, newV2, newV3, newV4, newV5, newV6, newV7);
    }
    
    /**
     * Rotates the cube around the X-axis by the specified angle.
     * 
     * This rotation pivots around the X-axis passing through the cube's center.
     * Follows right-hand rule: positive angles rotate counterclockwise when
     * looking along the positive X-axis toward the origin.
     * 
     * <p><strong>Algorithm:</strong></p>
     * <ol>
     *   <li>Translate cube to origin (center at 0,0,0)</li>
     *   <li>Apply rotation to each vertex using rotation matrix</li>
     *   <li>Translate back to original center position</li>
     * </ol>
     * 
     * <p><strong>Applications:</strong> Object orientation in 3D modeling, gimbal
     * rotations in flight simulators, animated transformations, user-controlled rotation.</p>
     * 
     * <p><strong>Important Note:</strong> After rotation, the cube may no longer be
     * axis-aligned, affecting efficiency of some spatial queries.</p>
     * 
     * <p><strong>Complexity:</strong> O(1) time (8 vertex rotations), O(1) space</p>
     * 
     * @param angleRadians the rotation angle in radians (not degrees)
     * @return a new Cube3D rotated around the X-axis
     */
    public Cube3D rotateX(double angleRadians) {
        if (Double.isNaN(angleRadians) || Double.isInfinite(angleRadians)) {
            logger.warning("Invalid rotation angle: " + angleRadians + ". Returning original cube.");
            return this;
        }
        
        Point3D center = getCenter();
        
        Point3D newV0 = translateToOrigin(v0, center).rotateX(angleRadians).translate(
            center.getX(), center.getY(), center.getZ());
        Point3D newV1 = translateToOrigin(v1, center).rotateX(angleRadians).translate(
            center.getX(), center.getY(), center.getZ());
        Point3D newV2 = translateToOrigin(v2, center).rotateX(angleRadians).translate(
            center.getX(), center.getY(), center.getZ());
        Point3D newV3 = translateToOrigin(v3, center).rotateX(angleRadians).translate(
            center.getX(), center.getY(), center.getZ());
        Point3D newV4 = translateToOrigin(v4, center).rotateX(angleRadians).translate(
            center.getX(), center.getY(), center.getZ());
        Point3D newV5 = translateToOrigin(v5, center).rotateX(angleRadians).translate(
            center.getX(), center.getY(), center.getZ());
        Point3D newV6 = translateToOrigin(v6, center).rotateX(angleRadians).translate(
            center.getX(), center.getY(), center.getZ());
        Point3D newV7 = translateToOrigin(v7, center).rotateX(angleRadians).translate(
            center.getX(), center.getY(), center.getZ());
        
        logger.info("Rotated cube around X-axis by " + Math.toDegrees(angleRadians) + " degrees");
        
        return new Cube3D(newV0, newV1, newV2, newV3, newV4, newV5, newV6, newV7);
    }
    
    /**
     * Rotates the cube around the Y-axis by the specified angle.
     * 
     * This rotation pivots around the Y-axis passing through the cube's center.
     * Follows right-hand rule: positive angles rotate counterclockwise when looking
     * along the positive Y-axis toward the origin.
     * 
     * <p>See {@link #rotateX(double)} for detailed explanation of the rotation
     * algorithm and applications.</p>
     * 
     * @param angleRadians the rotation angle in radians (not degrees)
     * @return a new Cube3D rotated around the Y-axis
     */
    public Cube3D rotateY(double angleRadians) {
        if (Double.isNaN(angleRadians) || Double.isInfinite(angleRadians)) {
            logger.warning("Invalid rotation angle: " + angleRadians + ". Returning original cube.");
            return this;
        }
        
        Point3D center = getCenter();
        
        Point3D newV0 = translateToOrigin(v0, center).rotateY(angleRadians).translate(
            center.getX(), center.getY(), center.getZ());
        Point3D newV1 = translateToOrigin(v1, center).rotateY(angleRadians).translate(
            center.getX(), center.getY(), center.getZ());
        Point3D newV2 = translateToOrigin(v2, center).rotateY(angleRadians).translate(
            center.getX(), center.getY(), center.getZ());
        Point3D newV3 = translateToOrigin(v3, center).rotateY(angleRadians).translate(
            center.getX(), center.getY(), center.getZ());
        Point3D newV4 = translateToOrigin(v4, center).rotateY(angleRadians).translate(
            center.getX(), center.getY(), center.getZ());
        Point3D newV5 = translateToOrigin(v5, center).rotateY(angleRadians).translate(
            center.getX(), center.getY(), center.getZ());
        Point3D newV6 = translateToOrigin(v6, center).rotateY(angleRadians).translate(
            center.getX(), center.getY(), center.getZ());
        Point3D newV7 = translateToOrigin(v7, center).rotateY(angleRadians).translate(
            center.getX(), center.getY(), center.getZ());
        
        logger.info("Rotated cube around Y-axis by " + Math.toDegrees(angleRadians) + " degrees");
        
        return new Cube3D(newV0, newV1, newV2, newV3, newV4, newV5, newV6, newV7);
    }
    
    /**
     * Rotates the cube around the Z-axis by the specified angle.
     * 
     * This rotation pivots around the Z-axis passing through the cube's center.
     * Follows right-hand rule: positive angles rotate counterclockwise when looking
     * along the positive Z-axis toward the origin.
     * 
     * <p>See {@link #rotateX(double)} for detailed explanation of the rotation
     * algorithm and applications.</p>
     * 
     * @param angleRadians the rotation angle in radians (not degrees)
     * @return a new Cube3D rotated around the Z-axis
     */
    public Cube3D rotateZ(double angleRadians) {
        if (Double.isNaN(angleRadians) || Double.isInfinite(angleRadians)) {
            logger.warning("Invalid rotation angle: " + angleRadians + ". Returning original cube.");
            return this;
        }
        
        Point3D center = getCenter();
        
        Point3D newV0 = translateToOrigin(v0, center).rotateZ(angleRadians).translate(
            center.getX(), center.getY(), center.getZ());
        Point3D newV1 = translateToOrigin(v1, center).rotateZ(angleRadians).translate(
            center.getX(), center.getY(), center.getZ());
        Point3D newV2 = translateToOrigin(v2, center).rotateZ(angleRadians).translate(
            center.getX(), center.getY(), center.getZ());
        Point3D newV3 = translateToOrigin(v3, center).rotateZ(angleRadians).translate(
            center.getX(), center.getY(), center.getZ());
        Point3D newV4 = translateToOrigin(v4, center).rotateZ(angleRadians).translate(
            center.getX(), center.getY(), center.getZ());
        Point3D newV5 = translateToOrigin(v5, center).rotateZ(angleRadians).translate(
            center.getX(), center.getY(), center.getZ());
        Point3D newV6 = translateToOrigin(v6, center).rotateZ(angleRadians).translate(
            center.getX(), center.getY(), center.getZ());
        Point3D newV7 = translateToOrigin(v7, center).rotateZ(angleRadians).translate(
            center.getX(), center.getY(), center.getZ());
        
        logger.info("Rotated cube around Z-axis by " + Math.toDegrees(angleRadians) + " degrees");
        
        return new Cube3D(newV0, newV1, newV2, newV3, newV4, newV5, newV6, newV7);
    }
    
    /**
     * Scales the cube by the specified factor around its center.
     * 
     * Scaling multiplies the distance of each vertex from the center by the scale
     * factor, making the cube larger (factor > 1) or smaller (0 < factor < 1).
     * 
     * <p><strong>Algorithm:</strong></p>
     * <ol>
     *   <li>Calculate cube center</li>
     *   <li>For each vertex: vector = (vertex - center) × scale</li>
     *   <li>New vertex = center + vector</li>
     * </ol>
     * 
     * <p><strong>Scaling Properties:</strong></p>
     * <ul>
     *   <li>Volume scales by factor³</li>
     *   <li>Surface area scales by factor²</li>
     *   <li>Edge lengths scale by factor</li>
     *   <li>Center position remains unchanged</li>
     * </ul>
     * 
     * <p><strong>Applications:</strong> LOD systems, object resizing in 3D editors,
     * zooming effects, bounding box expansion for collision margins.</p>
     * 
     * <p><strong>Complexity:</strong> O(1) time, O(1) space</p>
     * 
     * @param factor the scaling factor; values > 1 enlarge, 0 < values < 1 shrink
     * @return a new Cube3D scaled by the specified factor
     * @throws IllegalArgumentException if factor is non-positive
     */
    public Cube3D scale(double factor) {
        if (factor <= 0) {
            logger.severe("Attempted to scale cube by non-positive factor: " + factor);
            throw new IllegalArgumentException("Scale factor must be positive");
        }
        
        if (Double.isNaN(factor) || Double.isInfinite(factor)) {
            logger.warning("Invalid scale factor: " + factor + ". Returning original cube.");
            return this;
        }
        
        Point3D center = getCenter();
        
        Point3D newV0 = scalePointFromCenter(v0, center, factor);
        Point3D newV1 = scalePointFromCenter(v1, center, factor);
        Point3D newV2 = scalePointFromCenter(v2, center, factor);
        Point3D newV3 = scalePointFromCenter(v3, center, factor);
        Point3D newV4 = scalePointFromCenter(v4, center, factor);
        Point3D newV5 = scalePointFromCenter(v5, center, factor);
        Point3D newV6 = scalePointFromCenter(v6, center, factor);
        Point3D newV7 = scalePointFromCenter(v7, center, factor);
        
        logger.info("Scaled cube by factor " + factor);
        
        return new Cube3D(newV0, newV1, newV2, newV3, newV4, newV5, newV6, newV7);
    }
    
    /**
     * Tests whether a point is contained within this cube.
     * 
     * For an axis-aligned cube, a point is inside if its coordinates fall within
     * the ranges defined by the cube's extents on each axis. This is one of the most
     * fundamental spatial queries in 3D graphics and computational geometry.
     * 
     * <p><strong>Algorithm:</strong> Finds min/max coordinates along each axis, then
     * checks if the point's coordinates fall within all three ranges.</p>
     * 
     * <p><strong>Boundary Handling:</strong> Points exactly on the cube's surface are
     * considered inside (inclusive boundaries).</p>
     * 
     * <p><strong>Applications:</strong> Frustum culling (is object visible?), click/pick
     * testing in 3D editors, spatial partitioning and octree node assignment, collision
     * detection broad phase, region-based queries.</p>
     * 
     * <p><strong>Complexity:</strong> O(1) time for axis-aligned cubes.</p>
     * 
     * @param point the point to test; must not be null
     * @return true if the point is inside or on the boundary of the cube
     * @throws IllegalArgumentException if point is null
     */
    public boolean contains(Point3D point) {
        if (point == null) {
            logger.severe("Attempted to test containment of null point");
            throw new IllegalArgumentException("Cannot test containment of null point");
        }
        
        double minX = Math.min(v0.getX(), v6.getX());
        double maxX = Math.max(v0.getX(), v6.getX());
        double minY = Math.min(v0.getY(), v6.getY());
        double maxY = Math.max(v0.getY(), v6.getY());
        double minZ = Math.min(v0.getZ(), v6.getZ());
        double maxZ = Math.max(v0.getZ(), v6.getZ());
        
        boolean inside = point.getX() >= minX && point.getX() <= maxX &&
                        point.getY() >= minY && point.getY() <= maxY &&
                        point.getZ() >= minZ && point.getZ() <= maxZ;
        
        logger.info("Point " + point + " is " + (inside ? "inside" : "outside") + " cube");
        
        return inside;
    }
    
    /**
     * Tests whether this cube intersects with another cube.
     * 
     * Two axis-aligned cubes intersect if their projections onto all three axes overlap.
     * This is the Separating Axis Theorem (SAT) applied to AABBs - one of the most
     * efficient broad-phase collision detection algorithms.
     * 
     * <p><strong>Algorithm:</strong> For each axis (X, Y, Z), checks if the intervals
     * [min, max] of both cubes overlap. Cubes intersect only if they overlap on all axes.</p>
     * 
     * <p><strong>Interval Overlap Test:</strong> Two intervals [a, b] and [c, d] overlap
     * if: a ≤ d AND c ≤ b</p>
     * 
     * <p><strong>Applications:</strong> Broad-phase collision detection in physics engines,
     * frustum culling optimization, spatial index queries (octrees, k-d trees), region
     * intersection tests.</p>
     * 
     * <p><strong>Complexity:</strong> O(1) time - three axis comparisons regardless of
     * cube sizes or positions.</p>
     * 
     * @param other the other cube to test for intersection; must not be null
     * @return true if the cubes intersect (including touching at boundaries)
     * @throws IllegalArgumentException if other is null
     */
    public boolean intersects(Cube3D other) {
        if (other == null) {
            logger.severe("Attempted to test intersection with null cube");
            throw new IllegalArgumentException("Cannot test intersection with null cube");
        }
        
        double minX1 = Math.min(v0.getX(), v6.getX());
        double maxX1 = Math.max(v0.getX(), v6.getX());
        double minY1 = Math.min(v0.getY(), v6.getY());
        double maxY1 = Math.max(v0.getY(), v6.getY());
        double minZ1 = Math.min(v0.getZ(), v6.getZ());
        double maxZ1 = Math.max(v0.getZ(), v6.getZ());
        
        double minX2 = Math.min(other.v0.getX(), other.v6.getX());
        double maxX2 = Math.max(other.v0.getX(), other.v6.getX());
        double minY2 = Math.min(other.v0.getY(), other.v6.getY());
        double maxY2 = Math.max(other.v0.getY(), other.v6.getY());
        double minZ2 = Math.min(other.v0.getZ(), other.v6.getZ());
        double maxZ2 = Math.max(other.v0.getZ(), other.v6.getZ());
        
        boolean overlapX = minX1 <= maxX2 && minX2 <= maxX1;
        boolean overlapY = minY1 <= maxY2 && minY2 <= maxY1;
        boolean overlapZ = minZ1 <= maxZ2 && minZ2 <= maxZ1;
        
        boolean intersects = overlapX && overlapY && overlapZ;
        
        logger.info("Cubes " + (intersects ? "intersect" : "do not intersect"));
        
        return intersects;
    }
    
    /**
     * Returns an unmodifiable list of all vertices in the cube.
     * 
     * This method provides access to all eight vertices in their canonical ordering
     * (v0 through v7). The returned list is unmodifiable to preserve immutability.
     * 
     * <p><strong>Design Pattern:</strong> Collection Pattern with defensive copying
     * to prevent external modification while allowing traversal.</p>
     * 
     * <p><strong>Applications:</strong> Rendering vertex buffers for graphics APIs,
     * computing bounding volumes, mesh export to file formats, vertex-based spatial queries.</p>
     * 
     * @return an unmodifiable list of the eight vertices
     */
    public List<Point3D> getVertices() {
        List<Point3D> vertices = new ArrayList<>(8);
        vertices.add(v0);
        vertices.add(v1);
        vertices.add(v2);
        vertices.add(v3);
        vertices.add(v4);
        vertices.add(v5);
        vertices.add(v6);
        vertices.add(v7);
        
        return Collections.unmodifiableList(vertices);
    }
    
    /**
     * Returns an unmodifiable list of all 12 edges in the cube.
     * 
     * A cube has 12 edges connecting its 8 vertices: 4 on the front face, 4 on the
     * back face, and 4 connecting front to back. This method creates Line3D objects
     * for each edge following a consistent connectivity pattern.
     * 
     * <p><strong>Edge Connectivity:</strong></p>
     * <ul>
     *   <li>Front face (4 edges): v0-v1, v1-v2, v2-v3, v3-v0</li>
     *   <li>Back face (4 edges): v4-v5, v5-v6, v6-v7, v7-v4</li>
     *   <li>Connecting edges (4 edges): v0-v4, v1-v5, v2-v6, v3-v7</li>
     * </ul>
     * 
     * <p><strong>Applications:</strong> Wireframe rendering in 3D graphics, edge-based
     * collision detection, graph representations of mesh topology, CAD feature extraction.</p>
     * 
     * <p><strong>Data Structure Insight:</strong> This demonstrates the graph structure
     * of a cube - 8 vertices (nodes) connected by 12 edges with degree 3 at each vertex.
     * Understanding this topology is fundamental to mesh processing algorithms.</p>
     * 
     * @return an unmodifiable list of the twelve edges
     */
    public List<Line3D> getEdges() {
        List<Line3D> edges = new ArrayList<>(12);
        
        // Front face edges (Z-)
        edges.add(new Line3D(v0, v1));
        edges.add(new Line3D(v1, v2));
        edges.add(new Line3D(v2, v3));
        edges.add(new Line3D(v3, v0));
        
        // Back face edges (Z+)
        edges.add(new Line3D(v4, v5));
        edges.add(new Line3D(v5, v6));
        edges.add(new Line3D(v6, v7));
        edges.add(new Line3D(v7, v4));
        
        // Connecting edges (front to back)
        edges.add(new Line3D(v0, v4));
        edges.add(new Line3D(v1, v5));
        edges.add(new Line3D(v2, v6));
        edges.add(new Line3D(v3, v7));
        
        logger.info("Generated 12 edges for cube");
        
        return Collections.unmodifiableList(edges);
    }
    
    /**
     * Calculates the minimum distance from a point to the cube's surface.
     * 
     * This method computes the shortest distance from any point in 3D space to the
     * closest point on the cube's surface. For points inside the cube, the distance is zero.
     * 
     * <p><strong>Algorithm:</strong></p>
     * <ul>
     *   <li>If point is inside: distance = 0</li>
     *   <li>If point is outside: clamp the point to the cube's bounding box and
     *       compute distance to the clamped point</li>
     * </ul>
     * 
     * <p><strong>Applications:</strong> Proximity queries in physics engines, influence
     * radius calculations, distance field generation for rendering, nearest-object
     * queries in spatial databases.</p>
     * 
     * <p><strong>Complexity:</strong> O(1) time, O(1) space</p>
     * 
     * @param point the point to measure distance from; must not be null
     * @return the minimum distance to the cube's surface; 0 if point is inside
     * @throws IllegalArgumentException if point is null
     */
    public double distanceToPoint(Point3D point) {
        if (point == null) {
            logger.severe("Attempted to calculate distance to null point");
            throw new IllegalArgumentException("Cannot calculate distance to null point");
        }
        
        if (contains(point)) {
            logger.info("Point " + point + " is inside cube; distance = 0");
            return 0.0;
        }
        
        double minX = Math.min(v0.getX(), v6.getX());
        double maxX = Math.max(v0.getX(), v6.getX());
        double minY = Math.min(v0.getY(), v6.getY());
        double maxY = Math.max(v0.getY(), v6.getY());
        double minZ = Math.min(v0.getZ(), v6.getZ());
        double maxZ = Math.max(v0.getZ(), v6.getZ());
        
        double clampedX = Math.max(minX, Math.min(point.getX(), maxX));
        double clampedY = Math.max(minY, Math.min(point.getY(), maxY));
        double clampedZ = Math.max(minZ, Math.min(point.getZ(), maxZ));
        
        Point3D closestPoint = new Point3D(clampedX, clampedY, clampedZ);
        double distance = point.distanceTo(closestPoint);
        
        logger.info("Distance from point " + point + " to cube surface: " + distance);
        
        return distance;
    }
    
    // Helper methods
    
    private Point3D translateToOrigin(Point3D point, Point3D center) {
        return new Point3D(
            point.getX() - center.getX(),
            point.getY() - center.getY(),
            point.getZ() - center.getZ()
        );
    }
    
    private Point3D scalePointFromCenter(Point3D point, Point3D center, double factor) {
        double dx = (point.getX() - center.getX()) * factor;
        double dy = (point.getY() - center.getY()) * factor;
        double dz = (point.getZ() - center.getZ()) * factor;
        
        return new Point3D(
            center.getX() + dx,
            center.getY() + dy,
            center.getZ() + dz
        );
    }
    
    @Override
    public boolean equals(Object obj) {
        if (this == obj) return true;
        if (obj == null || getClass() != obj.getClass()) return false;
        
        Cube3D other = (Cube3D) obj;
        
        return v0.equals(other.v0) && v1.equals(other.v1) &&
               v2.equals(other.v2) && v3.equals(other.v3) &&
               v4.equals(other.v4) && v5.equals(other.v5) &&
               v6.equals(other.v6) && v7.equals(other.v7);
    }
    
    @Override
    public int hashCode() {
        int result = 17;
        result = 31 * result + v0.hashCode();
        result = 31 * result + v1.hashCode();
        result = 31 * result + v2.hashCode();
        result = 31 * result + v3.hashCode();
        result = 31 * result + v4.hashCode();
        result = 31 * result + v5.hashCode();
        result = 31 * result + v6.hashCode();
        result = 31 * result + v7.hashCode();
        return result;
    }
    
    @Override
    public String toString() {
        Point3D center = getCenter();
        double[] dims = getDimensions();
        return String.format("Cube3D[center=%s, dimensions=(%.3f, %.3f, %.3f)]",
                           center, dims[0], dims[1], dims[2]);
    }
}