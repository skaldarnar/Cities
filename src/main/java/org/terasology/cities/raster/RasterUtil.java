/*
 * Copyright 2015 MovingBlocks
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package org.terasology.cities.raster;

import org.terasology.math.TeraMath;
import org.terasology.math.geom.LineSegment;
import org.terasology.math.geom.Rect2i;
import org.terasology.math.geom.Vector2f;

/**
 * Converts model elements into blocks
 */
public abstract class RasterUtil {

//    /**
//     * @param shape the shape to fill
//     * @param hmBottom the bottom height map (inclusive)
//     * @param hmTop top height map (exclusive)
//     * @param type the block type
//     */
//    public static void fillShape(Shape shape, HeightMap hmBottom, HeightMap hmTop, BlockTypes type) {
//
//        if (!getAffectedArea().overlaps(shape.getBounds()) {
//            return;
//        }
//
//        Rect2i rc = getIntersectionArea(shape.getBounds());
//
//        for (int z = rc.minY(); z <= rc.maxY(); z++) {
//            for (int x = rc.minX(); x <= rc.maxX(); x++) {
//
//                if (shape.contains(x, z)) {
//                    int y1 = hmBottom.apply(x, z);
//                    int y2 = hmTop.apply(x, z);
//
//                    for (int y = y1; y < y2; y++) {
//                        setBlock(x, y, z, type);
//                    }
//                }
//            }
//        }
//    }

    /**
     * if (x2 < x1) nothing will be drawn.
     *
     * @param pen the pen use
     * @param x1  left x coord
     * @param x2  right x coord
     * @param z   the z coord
     */
    public static void drawLineX(Pen pen, int x1, int x2, int z) {
        Rect2i rc = pen.getTargetArea();

        if (z >= rc.minY() && z <= rc.maxY()) {
            int minX = Math.max(x1, rc.minX());
            int maxX = Math.min(x2, rc.maxX());
            for (int x = minX; x <= maxX; x++) {
                pen.draw(x, z);
            }
        }
    }

    /**
     * if (z2 < z1) nothing will be drawn.
     *
     * @param pen the pen use
     * @param z1  top z coord
     * @param z2  bottom z coord
     * @param x   the x coord
     */
    public static void drawLineZ(Pen pen, int x, int z1, int z2) {
        Rect2i rc = pen.getTargetArea();

        if (x >= rc.minX() && x <= rc.maxX()) {
            int minZ = Math.max(z1, rc.minY());
            int maxZ = Math.min(z2, rc.maxY());
            for (int z = minZ; z <= maxZ; z++) {
                pen.draw(x, z);
            }
        }
    }

    /**
     * @param rect the area to fill
     * @param pen  the pen to use for the rasterization of the rectangle
     */
    public static void fillRect(Pen pen, Rect2i rect) {
        Rect2i rc = pen.getTargetArea().intersect(rect);

        if (rc.isEmpty()) {
            return;
        }

        for (int z = rc.minY(); z <= rc.maxY(); z++) {
            for (int x = rc.minX(); x <= rc.maxX(); x++) {
                pen.draw(x, z);
            }
        }
    }

    /**
     * @param pen the pen to use
     * @param rc  the rectangle to draw
     */
    public static void drawRect(Pen pen, Rect2i rc) {

        // walls along x-axis
        drawLineX(pen, rc.minX(), rc.maxX(), rc.minY());
        drawLineX(pen, rc.minX(), rc.maxX(), rc.maxY());

        // walls along z-axis
        drawLineZ(pen, rc.minX(), rc.minY() + 1, rc.maxY() - 1); // no need to draw corners again
        drawLineZ(pen, rc.maxX(), rc.minY() + 1, rc.maxY() - 1); //  -> inset by one on both ends

    }

    /**
     * Draws a line.<br>
     * See Wikipedia: Bresenham's line algorithm, chapter Simplification
     *
     * @param pen  the pen to use
     * @param line the line to draw
     */
    public static void drawLine(Pen pen, LineSegment line) {

        Rect2i outerBox = pen.getTargetArea();

        Vector2f p0 = new Vector2f();
        Vector2f p1 = new Vector2f();
        if (line.getClipped(outerBox, p0, p1)) {
            int cx1 = TeraMath.floorToInt(p0.getX());
            int cy1 = TeraMath.floorToInt(p0.getY());
            int cx2 = TeraMath.floorToInt(p1.getX());
            int cy2 = TeraMath.floorToInt(p1.getY());
            drawClippedLine(pen, cx1, cy1, cx2, cy2);
        }
    }

    private static void drawClippedLine(Pen pen, int x1, int z1, int x2, int z2) {

        int dx = Math.abs(x2 - x1);
        int dy = Math.abs(z2 - z1);

        int sx = (x1 < x2) ? 1 : -1;
        int sy = (z1 < z2) ? 1 : -1;

        int err = dx - dy;

        int x = x1;
        int z = z1;

        while (true) {
            pen.draw(x, z);

            if (x == x2 && z == z2) {
                break;
            }

            int e2 = 2 * err;

            if (e2 > -dy) {
                err = err - dy;
                x += sx;
            }
            // if going along diagonals is not ok use " .. } else if (e2.. " instead

            if (e2 < dx) {
                err = err + dx;
                z += sy;
            }
        }
    }

    public enum LineOverlap {
        NONE,
        MAJOR,
        MINOR,
        BOTH

    }

    public enum ThicknessMode {
        MIDDLE,
        CLOCKWISE,
        COUNTERCLOCKWISE

    }

    /**
     * Modified Bresenham with optional overlap (esp. for drawThickLine())
     * Overlap draws additional pixel when changing minor direction - for standard bresenham overlap = LINE_OVERLAP_NONE (0)
     * <p>
     * Sample line:
     * <p>
     * 00+
     * -0000+
     * -0000+
     * -00
     * <p>
     * 0 pixels are drawn for normal line without any overlap LINE_OVERLAP_NONE
     * + pixels are drawn if LINE_OVERLAP_MAJOR
     * - pixels are drawn if LINE_OVERLAP_MINOR
     *
     * @param pen
     * @param x1      start x coordinate
     * @param y1      start z coordinate
     * @param x2      end x coordinate
     * @param y2      end z coordinate
     * @param overlap the overlap mode
     */
    public static void drawLineOverlap(Pen pen, int x1, int y1, int x2, int y2, LineOverlap overlap) {

        // horizontal or vertical line -> fillRect() is faster
        if (x1 == x2 || y1 == y2) {
            fillRect(pen, Rect2i.createFromMinAndMax(x1, y1, x2, y2));
            return;
        }

        int dx = Math.abs(x2 - x1);
        int dy = Math.abs(y2 - y1);

        int sx = (x1 < x2) ? 1 : -1;
        int sy = (y1 < y2) ? 1 : -1;

        int error;

        // draw first pixel
        pen.draw(x1, y1);
        if (dx > dy) {
            // start value represents a half step in Y direction
            error = (dy * 2) - dx;
            while (x1 != x2) {
                //step in main direction
                x1 += sx;
                if (error >= 0) {
                    if (overlap == LineOverlap.MAJOR || overlap == LineOverlap.BOTH) {
                        // draw pixel in main direction before changing
                        pen.draw(x1, y1);
                    }
                    // change Y
                    y1 += sy;
                    if (overlap == LineOverlap.MINOR || overlap == LineOverlap.BOTH) {
                        // draw pixel in minor direction before changing
                        pen.draw(x1 - sx, y1);
                    }
                    error -= (dx * 2);
                }
                error -= (dy * 2);
                pen.draw(x1, y1);
            }
        } else {
            error = (dx * 2) - dy;
            while (y1 != y2) {
                //step in main direction
                y1 += sy;
                if (error >= 0) {
                    if (overlap == LineOverlap.MAJOR || overlap == LineOverlap.BOTH) {
                        // draw pixel in main direction before changing
                        pen.draw(x1, y1);
                    }
                    // change X
                    x1 += sx;
                    if (overlap == LineOverlap.MINOR || overlap == LineOverlap.BOTH) {
                        // draw pixel in minor direction before changing
                        pen.draw(x1, y1 - sy);
                    }
                    error -= (dy * 2);
                }
                error -= (dx * 2);
                pen.draw(x1, y1);
            }
        }
    }

    /**
     * Bresenham with thickness.
     * <p>
     * No pixel is missed, and every pixel is just drawn once.
     *
     * @param pen
     * @param x1
     * @param y1
     * @param x2
     * @param y2
     * @param thickness
     * @param mode
     */
    public static void drawThickLine(Pen pen, int x1, int y1, int x2, int y2, int thickness, ThicknessMode mode) {
        if (thickness <= 1) {
            drawLineOverlap(pen, x1, y1, x2, y2, LineOverlap.NONE);
        }

        // TODO: clip to display size/visible area

        int dx = TeraMath.fastAbs(x2 - x1);
        int dy = TeraMath.fastAbs(y2 - y1);
        int sx = (x1 < x2) ? 1 : -1;
        int sy = (y1 < y2) ? 1 : -1;

        int error;

        boolean swap = (sx == sy);
        LineOverlap overlap;

        // adjust for right direction from of thickness from line origin
        int drawStartAdjustCount = thickness / 2;
        switch (mode) {
            case COUNTERCLOCKWISE:
                drawStartAdjustCount = thickness - 1;
                break;
            case CLOCKWISE:
                drawStartAdjustCount = 0;
                break;
        }

        // in which octant are we in now
        if (dx >= dy) {
            if (swap) {
                drawStartAdjustCount = (thickness - 1) - drawStartAdjustCount;
                sy = -sy;
            } else {
                sx = -sx;
            }

            // adjust draw start point
            error = (dy * 2) - dx;
            for (int i = drawStartAdjustCount; i > 0; i--) {
                // change X (main direction here)
                x1 -= sx;
                x2 -= sx;
                if (error >= 0) {
                    // change Y
                    y1 -= sy;
                    y2 -= sy;
                    error -= dx * 2;
                }
                error += dy * 2;
            }

            // draw start line
            drawClippedLine(pen, x1, y1, x2, y2);

            // draw `thickness` lines
            error = (dy * 2) - dx;
            for (int i = thickness; i > 1; i--) {
                // change X (main direction here)
                x1 += sx;
                x2 += sx;
                overlap = LineOverlap.NONE;
                if (error >= 0) {
                    // change Y
                    y1 += sy;
                    y2 += sy;
                    error -= dx * 2;
                    /*
                 * change in minor direction reverse to line (main) direction
				 * because of choosing the right (counter)clockwise draw vector
				 * use LINE_OVERLAP_MAJOR to fill all pixel
				 *
				 * EXAMPLE:
				 * 1,2 = Pixel of first lines
				 * 3 = Pixel of third line in normal line mode
				 * - = Pixel which will additionally be drawn in LINE_OVERLAP_MAJOR mode
				 *           33
				 *       3333-22
				 *   3333-222211
				 * 33-22221111
				 *  221111                     /\
				 *  11                          Main direction of draw vector
				 *  -> Line main direction
				 *  <- Minor direction of counterclockwise draw vector
				 */
                    overlap = LineOverlap.MAJOR;
                }
                error += dy * 2;
                drawLineOverlap(pen, x1, y1, x2, y2, overlap);
            }
        } else {
            // the other octant
            if (swap) {
                sx = -sx;
            } else {
                drawStartAdjustCount = (thickness - 1) - drawStartAdjustCount;
                sy = -sy;
            }
            // adjust draw start point
            error = dx * 2 - dy;
            for (int i = drawStartAdjustCount; i > 0; i--) {
                y1 -= sy;
                y2 -= sy;
                if (error >= 0) {
                    x1 -= sx;
                    x2 -= sx;
                    error -= dy * 2;
                }
                error += dx * 2;
            }
            //draw start line
            drawClippedLine(pen, x1, y1, x2, y2);
            error = dx * 2 - dy;
            for (int i = thickness; i > 1; i--) {
                y1 += sy;
                y2 += sy;
                overlap = LineOverlap.NONE;
                if (error >= 0) {
                    x1 += sx;
                    x2 += sx;
                    error -= dy * 2;
                    overlap = LineOverlap.MAJOR;
                }
                error += dx * 2;
                drawLineOverlap(pen, x1, y1, x2, y2, overlap);
            }
        }
    }


    public static void drawThickLineSimple(Pen pen, LineSegment line, float thickness) {
        Rect2i outerBox = pen.getTargetArea();

        Vector2f p0 = new Vector2f();
        Vector2f p1 = new Vector2f();
        if (line.getClipped(outerBox, p0, p1)) {
            int cx1 = TeraMath.floorToInt(p0.getX());
            int cy1 = TeraMath.floorToInt(p0.getY());
            int cx2 = TeraMath.floorToInt(p1.getX());
            int cy2 = TeraMath.floorToInt(p1.getY());
            drawClippedThickLineSimple(pen, cx1, cy1, cx2, cy2, thickness);
        }
    }

    public static void drawClippedThickLineSimple(Pen pen, int x0, int y0, int x1, int y1, float thickness) {
        int dx = TeraMath.fastAbs(x1 - x0);
        int dy = TeraMath.fastAbs(y1 - y0);
        int sx = x0 < x1 ? 1 : -1;
        int sy = y0 < y1 ? 1 : -1;

        int err = dx - dy;                        /* error value e_xy */
        int e2, x2, y2;
        float alpha;

        float ed = dx + dy == 0 ? 1 : TeraMath.sqrt((float) dx * dx + (float) dy * dy);

        for (thickness = (thickness + 1) / 2; ; ) {                                           /* pixel loop */
            alpha = Math.max(0, TeraMath.fastAbs(err - dx + dy) / ed - thickness + 1);
            if (alpha > 0.2) {
                pen.draw(x0, y0);
            }
            pen.draw(x0, y0);
            e2 = err;
            x2 = x0;
            if (2 * e2 >= -dx) {                                           /* x step */
                for (e2 += dy, y2 = y0; e2 < ed * thickness && (y1 != y2 || dx > dy); e2 += dx) {
                    alpha = Math.max(0, TeraMath.fastAbs(e2) / ed - thickness + 1);
                    if (alpha > 0.2) {
                        pen.draw(x0, y2);
                    }
                    pen.draw(x0, y2);
                    y2 += sy;
                }
                if (x0 == x1) break;
                e2 = err;
                err -= dy;
                x0 += sx;
            }
            if (2 * e2 <= dy) {                                            /* y step */
                for (e2 = dx - e2; e2 < ed * thickness && (x1 != x2 || dx < dy); e2 += dy) {
                    alpha = Math.max(0, TeraMath.fastAbs(e2) / ed - thickness + 1);
                    if (alpha > 0.2) {
                        pen.draw(x2, y0);
                    }
                    pen.draw(x2, y0);
                    x2 += sx;
                }
                if (y0 == y1) break;
                err += dx;
                y0 += sy;
            }
        }
    }

    /**
     * Draws a circle based on Horn's algorithm (see B. K. P. Horn: Circle Generators for Display Devices.
     * Computer Graphics and Image Processing 5, 2 - June 1976)
     *
     * @param cx         the center x
     * @param cy         the center y
     * @param rad        the radius
     * @param checkedPen the receiving instance. Must be checked, because iterator could draw outside.
     */
    public static void drawCircle(CheckedPen checkedPen, int cx, int cy, int rad) {
        int d = -rad;
        int x = rad;
        int y = 0;
        while (y <= x) {
            checkedPen.draw(cx + x, cy + y);
            checkedPen.draw(cx - x, cy + y);
            checkedPen.draw(cx - x, cy - y);
            checkedPen.draw(cx + x, cy - y);

            checkedPen.draw(cx + y, cy + x);
            checkedPen.draw(cx - y, cy + x);
            checkedPen.draw(cx - y, cy - x);
            checkedPen.draw(cx + y, cy - x);

            d = d + 2 * y + 1;
            y = y + 1;
            if (d > 0) {
                d = d - 2 * x + 2;
                x = x - 1;
            }
        }
    }

    /**
     * @param cx  the center x
     * @param cy  the center y
     * @param rad the radius
     * @param pen the pen to draw
     */
    public static void fillCircle(CheckedPen pen, int cx, int cy, int rad) {
        for (int y = 0; y <= rad; y++) {
            for (int x = 0; x * x + y * y <= (rad + 0.5) * (rad + 0.5); x++) {
                pen.draw(cx + x, cy + y);
                pen.draw(cx - x, cy + y);
                pen.draw(cx - x, cy - y);
                pen.draw(cx + x, cy - y);
            }
        }
    }
}
