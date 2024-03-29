-- | Double-precision BSpline smoothing using hmatrix

module Math.Spline.Smoothing
    ( ordinary, unweighted, weighted, covWeighted, basicSmooth
    -- * Functions for convenience
    , smoothCubic, smoothQuadratic, iSpecial
    ) where

import Control.Arrow
import Control.Monad
import Data.List (transpose, zip4)

import Data.Packed.Matrix
import Data.Packed.Vector
import qualified Data.Vector.Storable as SV
import qualified Data.Vector as BV
import Math.Polynomial hiding (x)
import Math.Spline.BSpline
import Math.Spline.BSpline.Reference
import Math.Spline.Class
import Math.Spline.ISpline
import Math.Spline.Knots hiding (fromList, knots)
import Numeric.Container
import Numeric.LinearAlgebra.Algorithms

import Debug.Trace

ordinary :: Int -> [Double] -> [(Double, Double)] -> BSpline SV.Vector Double
ordinary degree knotList dataPts = bSpline knots solution
  where
    splineBases = basisFunctions knots !! degree
    knots = mkKnots knotList
    (dataXs, dataYs) = unzip dataPts
    dataYC = fromList dataYs
    coeffs = fromColumns $ map (\f -> fromList $ map f dataXs) splineBases
    solution = coeffs <\> dataYC

-- | Same as ordinary.
unweighted :: Int -> [Double] -> [(Double, Double)] -> BSpline SV.Vector Double
unweighted = ordinary

-- | Weighted smoothing spline, with weights supplied per point and with no
-- covariance. The weight matrix is simply the diagonal matrix of weights
weighted :: Int -> [Double]
         -> [(Double, Double, Double)]
         -- ^ a list of (abscissa, ordinate, weight) tuples
         -> BSpline SV.Vector Double
weighted degree knotList dataPts = bSpline knots solution
  where
    splineBases = basisFunctions knots !! degree
    knots = mkKnots knotList
    (dataXs, dataYs, weights) = unzip3 dataPts
    weightMatrix = diag $ fromList weights
    dataYC = fromList dataYs
    coeffs = fromColumns $ map (\f -> fromList $ map f dataXs) splineBases
    solution = (trans coeffs <> weightMatrix <> coeffs) <\>
               (trans coeffs <> weightMatrix <> dataYC)

-- | Weighted smoothing spline with weights and covariance supplied as a
-- weight matrix
covWeighted :: Int -> [Double] -> [(Double, Double)] -> Matrix Double
            -> BSpline SV.Vector Double
covWeighted degree knotList dataPts weightMatrix = bSpline knots solution
  where
    splineBases = basisFunctions knots !! degree
    knots = mkKnots knotList
    (dataXs, dataYs) = unzip dataPts
    dataYC = fromList dataYs
    coeffs = fromColumns $ map (\f -> fromList $ map f dataXs) splineBases
    solution = (trans coeffs <> weightMatrix <> coeffs) <\>
               (trans coeffs <> weightMatrix <> dataYC)

-- | Smoothing spline based on the (degree - 1)st derivative of the spline.
basicSmooth :: Double -- ^ Smoothness coefficient -- A non-negative number
            -> Int -> [Double] -> [(Double, Double)] -> Matrix Double
            -> BSpline SV.Vector Double
basicSmooth smoothness degree knotList dataPts weightMatrix =
    bSpline knots solution
  where
    splineBases = basisFunctions knots !! degree
    derivBases = transpose $ map (map ((!! (degree - 1)) . iterate polyDeriv)) $
                 map (!! degree) (basisPolynomials knots)
    dKnots = distinctKnots knots
    integrands = [ zip4 dKnots (tail dKnots) x y
                   | x <- derivBases, y <- derivBases]
    integrals = map (sum . map f) integrands
      where
        f (start, end, poly1, poly2) = p end - p start
          where
            p = evalPoly . polyIntegral $ multPoly poly1 poly2
    integralMatrix = ((length derivBases) >< (length derivBases)) integrals
    knots = mkKnots knotList
    (dataXs, dataYs) = unzip dataPts
    dataYC = fromList dataYs
    coeffs = fromColumns $ map (\f -> fromList $ map f dataXs) splineBases
    solution = head . toColumns $ cholSH ((trans coeffs <> weightMatrix <> coeffs)
                          `add` scale smoothness integralMatrix)
               `cholSolve` asColumn (trans coeffs <> weightMatrix <> dataYC)

smoothCubic :: Double -- ^ Smoothness coefficient -- A non-negative number
            -> [(Double, Double, Double)] -> BSpline SV.Vector Double
smoothCubic smoothness dataPts =
    basicSmooth smoothness 3 ks (zip xs ys) (diag $ fromList zs)
  where
    (xs,ys,zs) = unzip3 dataPts
    ks = replicate 3 (head xs) ++ xs
         ++ replicate 3 (2 * last xs - (last (init xs)))

smoothQuadratic :: Double -- ^ Smoothness coefficient -- A non-negative number
                -> [(Double, Double, Double)] -> BSpline SV.Vector Double
smoothQuadratic smoothness dataPts =
    basicSmooth smoothness 2 ks (zip xs ys) (diag $ fromList zs)
  where
    (xs,ys,zs) = unzip3 dataPts
    ks = replicate 2 (head xs) ++ xs
         ++ replicate 2 (2 * last xs - (last (init xs)))

-- | Like ordinary least squares, but where the derivative is constrained to be
-- an iSpline where the sum of control points is 1
--
-- Todo: Think of a more general interface to allow such things
iSpecial :: Int -> [Double] -> [(Double, Double)] -> BSpline SV.Vector Double
iSpecial degree knotList dataPts = bSpline knots normalSol
  where
    splineBases = basisFunctions knots !! degree
    knots = mkKnots knotList
    (dataXs, dataYs) = unzip dataPts
    dataYC = fromList dataYs
    coeffs = fromColumns $ map (\f -> fromList $ map f dataXs) splineBases
    numCtrlPts = numKnots knots - degree - 1
    (iSums, iSum) = (init &&& last) $ map f [0..(numCtrlPts - 1)]
    normalizeM = fromLists . (++ [map ((/ iSum) . negate) iSums]) . toLists $ diag (SV.replicate (numCtrlPts - 1) 1)
    normalizeDataYC = zipVectorWith (-) dataYC (coeffs <> last1)
    last1 = SV.replicate numCtrlPts 0 SV.// [(numCtrlPts - 1, 1/iSum)]
    f i = BV.sum . controlPoints . toISpline . differentiateBSpline $ bSpline knots ctl
      where
        ctl = SV.replicate numCtrlPts 0 SV.// [(i, 1)]
    sol = (coeffs <> normalizeM) <\> normalizeDataYC
    normalSol = zipVectorWith (+) (normalizeM <> sol) last1
