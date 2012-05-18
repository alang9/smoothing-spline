module Graphics.Rendering.Chart.Spline
    ( plotBSpline
    ) where

import qualified Data.Map as M

import Data.Accessor
import Data.Colour
import Data.Colour.Names
import Graphics.Rendering.Chart
import Math.Spline
import Math.Spline.BSpline
import Math.Spline.Knots hiding (fromList)

plotBSpline :: BSpline Double -> [(Double, Double, Double)] -> Layout1 Double Double
plotBSpline spline weightedDataPts = layout
  where
    knotVec = knotVector spline
    knotMap = toMap knotVec
    firstKnot = fst $ M.findMin knotMap
    lastKnot = fst $ M.findMax knotMap
    spline1 = plot_lines_values ^= [[(x, evalBSpline spline x) | x <- enumFromThenTo firstKnot (firstKnot + 0.01) lastKnot]]
              $ plot_lines_style .> line_color ^= opaque blue
              $ defaultPlotLines
    spline2 = plot_points_values ^= [(x, evalBSpline spline x) | x <- knots knotVec]
              $ plot_points_style ^= hollowCircles 4 1 (opaque blue)
              $ defaultPlotPoints
    spline3 = plot_lines_values ^= [[(x, evalBSpline (differentiateBSpline spline) x) | x <- enumFromThenTo firstKnot (firstKnot + 0.01) lastKnot]]
              $ plot_lines_style .> line_color ^= opaque purple
              $ defaultPlotLines
    spline4 = plot_lines_values ^= [[(x, evalBSpline (differentiateBSpline $ differentiateBSpline spline) x) | x <- enumFromThenTo firstKnot (firstKnot + 0.01) lastKnot]]
              $ plot_lines_style .> line_color ^= opaque pink
              $ defaultPlotLines
    dataPts1 = area_spots_values ^= weightedDataPts
              $ area_spots_fillcolour ^= flip withOpacity 0.2 green
              $ defaultAreaSpots
    dataPts2 = plot_points_values ^= map (\(x,y,_) -> (x,y)) weightedDataPts
              $ plot_points_style ^= filledCircles 2 (opaque green)
              $ defaultPlotPoints
    dataPts3 = plot_lines_values ^= map (\(x,y,_) -> [(x, y), (x, evalBSpline spline x)]) weightedDataPts
               $ plot_lines_style .> line_color ^= opaque red
               $ defaultPlotLines
    layout = layout1_plots ^= [ Left (toPlot spline1)
                              , Left (toPlot spline2)
                              , Left (toPlot dataPts1)
                              , Left (toPlot dataPts2)
                              , Left (toPlot dataPts3)
                              , Left (toPlot spline3)
                              , Left (toPlot spline4)
                              ] $ defaultLayout1
