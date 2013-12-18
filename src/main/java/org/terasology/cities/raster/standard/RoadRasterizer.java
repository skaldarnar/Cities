/*
 * Copyright 2013 MovingBlocks
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package org.terasology.cities.raster.standard;

import java.awt.Shape;

import org.terasology.cities.BlockTypes;
import org.terasology.cities.raster.Brush;
import org.terasology.cities.raster.Rasterizer;
import org.terasology.cities.raster.TerrainInfo;

/**
 * Draws road shapes on the terrain surface
 * @author Martin Steiger
 */
public class RoadRasterizer implements Rasterizer<Shape> {

    @Override
    public void raster(Brush brush, TerrainInfo ti, Shape element) {
        brush.fillShape(element, ti.getHeightMap(), 1, BlockTypes.ROAD_SURFACE);
    }

}
