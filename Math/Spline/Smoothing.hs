-- | Double-precision BSpline smoothing using hmatrix

module Math.Spline.Smoothing
    ( ordinary, unweighted, weighted, covWeighted, basicSmooth, smoothCubic
    ) where

import Data.List (transpose, zip4)

import Data.Packed.Matrix
import Data.Packed.Vector
import Data.Vector.Generic (convert)
import Data.Vector.Storable ()
import Math.Polynomial hiding (x)
import Math.Spline.BSpline
import Math.Spline.BSpline.Reference
import Math.Spline.Knots hiding (fromList, knots)
import Numeric.Container

ordinary :: Int -> [Double] -> [(Double, Double)] -> BSpline Double
ordinary degree knotList dataPts = bSpline knots solution
  where
    splineBases = basisFunctions knots !! degree
    knots = mkKnots knotList
    (dataXs, dataYs) = unzip dataPts
    dataYC = fromList dataYs
    coeffs = fromColumns $ map (\f -> fromList $ map f dataXs) splineBases
    solution = convert $ coeffs <\> dataYC

-- | Same as ordinary.
unweighted :: Int -> [Double] -> [(Double, Double)] -> BSpline Double
unweighted = ordinary

-- | Weighted smoothing spline, with weights supplied per point and with no
-- covariance. The weight matrix is simply the diagonal matrix of weights
weighted :: Int -> [Double]
         -> [(Double, Double, Double)]
         -- ^ a list of (abscissa, ordinate, weight) tuples
         -> BSpline Double
weighted degree knotList dataPts = bSpline knots solution
  where
    splineBases = basisFunctions knots !! degree
    knots = mkKnots knotList
    (dataXs, dataYs, weights) = unzip3 dataPts
    weightMatrix = diag $ fromList weights
    dataYC = fromList dataYs
    coeffs = fromColumns $ map (\f -> fromList $ map f dataXs) splineBases
    solution = convert $ (trans coeffs <> weightMatrix <> coeffs) <\>
               (trans coeffs <> weightMatrix <> dataYC)

-- | Weighted smoothing spline with weights and covariance supplied as a
-- weight matrix
covWeighted :: Int -> [Double] -> [(Double, Double)] -> Matrix Double
            -> BSpline Double
covWeighted degree knotList dataPts weightMatrix = bSpline knots solution
  where
    splineBases = basisFunctions knots !! degree
    knots = mkKnots knotList
    (dataXs, dataYs) = unzip dataPts
    dataYC = fromList dataYs
    coeffs = fromColumns $ map (\f -> fromList $ map f dataXs) splineBases
    solution = convert $ (trans coeffs <> weightMatrix <> coeffs) <\>
               (trans coeffs <> weightMatrix <> dataYC)

-- | Smoothing spline based on the (degree - 1)st derivative of the spline.
basicSmooth :: Double -> Int -> [Double] -> [(Double, Double)] -> Matrix Double
            -> BSpline Double
basicSmooth smoothness degree knotList dataPts weightMatrix =
    bSpline knots solution
  where
    splineBases = basisFunctions knots !! degree
    derivBases = transpose $ map (map ((!! (degree - 1)) . iterate polyDeriv)) $
                 map (!! degree) (basisPolynomials knots)
    dKnots = distinctKnots knots
    integrands = [ zip4 dKnots (tail dKnots) x y
                   | x <- derivBases, y <- derivBases]
    integrals = flip map integrands $ \l -> sum $ map f l
      where
        f (start, end, poly1, poly2) = p end - p start
          where
            p = evalPoly . polyIntegral $ multPoly poly1 poly2
    integralMatrix = ((length derivBases) >< (length derivBases)) integrals
    knots = mkKnots knotList
    (dataXs, dataYs) = unzip dataPts
    dataYC = fromList dataYs
    coeffs = fromColumns $ map (\f -> fromList $ map f dataXs) splineBases
    solution = convert $ ((trans coeffs <> weightMatrix <> coeffs)
                          `add` scale smoothness integralMatrix)
               <\> (trans coeffs <> weightMatrix <> dataYC)

smoothCubic :: Double -> [(Double, Double, Double)] -> BSpline Double
smoothCubic smoothness dataPts =
    basicSmooth smoothness 3 ks (zip xs ys) (diag $ fromList zs)
  where
    (xs,ys,zs) = unzip3 dataPts
    ks = replicate 3 (head xs) ++ xs
         ++ replicate 3 (2 * last xs - (last (init xs)))