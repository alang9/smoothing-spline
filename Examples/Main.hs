{-# OPTIONS_GHC -Wall #-}

module Main where

import Graphics.Rendering.Chart
import Graphics.Rendering.Chart.Gtk
import Math.Spline.BSpline
import Numeric.Container

import Math.Spline.Smoothing
import Graphics.Rendering.Chart.Spline

import Data.Accessor
import Data.Colour.Names
import Data.Colour
import qualified Math.Spline.Knots as K
import Math.Spline
import qualified Data.Vector as V
import qualified Data.Vector.Storable as VS

weightedDataPts1 :: [(Double, Double, Double)]
weightedDataPts1 = [(x, if x == 0 then 0 else sin (recip x), exp x) | x <- [-5,-4.75..5]]

spline1 :: BSpline VS.Vector Double
spline1 = basicSmooth 1e-3 3 ([-5,-5,-5] ++ map (\(x,_,_) -> x) weightedDataPts1 ++ [6,6,6])
          (map (\(x,y,_) -> (x,y)) weightedDataPts1)
          (diag . fromList $ map (\(_,_,z) -> z) weightedDataPts1)

spline2 :: BSpline VS.Vector Double
spline2 = weighted 3 ([-5,-5,-5] ++ map (\(x,_,_) -> x) weightedDataPts1 ++ [6,7,8])
          weightedDataPts1

spline3 :: BSpline VS.Vector Double
spline3 = smoothCubic 1e-3 weightedDataPts1

ytm :: Double
ytm = 1/12



f x = (evalNaturalBSpline spline3 x - evalNaturalBSpline (differentiateBSpline spline3) x / sqrt ytm
      + evalNaturalBSpline (differentiateBSpline . differentiateBSpline $ spline3) x / sqrt ytm ^ 2
            - evalNaturalBSpline (differentiateBSpline . differentiateBSpline . differentiateBSpline $ spline3) x / sqrt ytm ^ 3) * k + evalNaturalBSpline pc x
  where
    k = exp (sqrt ytm * x)

h x = evalNaturalBSpline (differentiateBSpline . differentiateBSpline . differentiateBSpline $ spline3) x * k / sqrt ytm ^ 3
  where
    k = exp (sqrt ytm * x)

pc = bSpline (knotVector pc') (V.postscanl' (\s (x, c, c') -> s + (c - c') * exp (sqrt ytm * x) / sqrt ytm ^ 3) 0 $ V.zip3 (K.knotsVector $ knotVector pc') (controlPoints pc') (0 `V.cons` controlPoints pc'))
  where
    pc' = differentiateBSpline . differentiateBSpline . differentiateBSpline $ spline3

g x = evalNaturalBSpline spline3 x * k
  where
    k = exp (sqrt ytm * x)

main :: IO ()
main = renderableToWindow (toRenderable (plotFn f)) 640 480

plotFn :: (Double -> Double) -> Layout1 Double Double
plotFn fun = layout
  where
    spline1 = plot_lines_values ^= [[(x, fun x / 10) | x <- [-5,-4.99..5]]]
              $ plot_lines_style .> line_color ^= opaque blue
              $ defaultPlotLines
    spline2 = plot_lines_values ^= [[(x, g x / 10) | x <- [-5,-4.99..5]]]
              $ plot_lines_style .> line_color ^= opaque black
              $ defaultPlotLines
    spline3 = plot_lines_values ^= [[(x, h x / 100) | x <- [-5,-4.99..5]]]
              $ plot_lines_style .> line_color ^= opaque orange
              $ defaultPlotLines
    spline4 = plot_lines_values ^= [[(x, evalNaturalBSpline pc x / 10000) | x <- [-5,-4.99..5]]]
              $ plot_lines_style .> line_color ^= opaque green
              $ defaultPlotLines
    layout = layout1_plots ^= [ Left (toPlot spline1)
                              , Left (toPlot spline2)
                     --         , Left (toPlot spline3)
                            --  , Left (toPlot spline4)
                              ] $ defaultLayout1
