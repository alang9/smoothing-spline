module Main where

import Graphics.Rendering.Chart
import Graphics.Rendering.Chart.Gtk
import Math.Spline.BSpline
import Numeric.Container

import Math.Spline.Smoothing
import Graphics.Rendering.Chart.Spline

weightedDataPts1 :: [(Double, Double, Double)]
weightedDataPts1 = [(x, sin (-x^2), exp x) | x <- [-5,-4.75..5]]

spline1 :: BSpline Double
spline1 = basicSmooth 1e-3 3 ([-5,-5,-5] ++ map (\(x,_,_) -> x) weightedDataPts1 ++ [6,6,6])
          (map (\(x,y,_) -> (x,y)) weightedDataPts1)
          (diag . fromList $ map (\(_,_,z) -> z) weightedDataPts1)

spline2 :: BSpline Double
spline2 = weighted 3 ([-5,-5,-5] ++ map (\(x,_,_) -> x) weightedDataPts1 ++ [6,7,8])
          weightedDataPts1

spline3 :: BSpline Double
spline3 = smoothCubic 1e-3 weightedDataPts1

main :: IO ()
main = renderableToWindow (toRenderable (plotBSpline spline3 weightedDataPts1)) 640 480
