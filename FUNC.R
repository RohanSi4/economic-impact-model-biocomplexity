###############################################################################
# Production: Calculates production output based on available stock,
# efficiency factors, and order constraints.
#
# Parameters:
#   StockMat  - Matrix representing the available stock for each sector-region.
#   E_CZ      - Efficiency (conversion) factors for the stock.
#   IOV_t     - Matrix representing inflows/outflows at time t.
#   E_VA      - Efficiency factor for the inflows/outflows.
#   OrderMat  - Matrix of orders for each sector-region.
#   IOX_0     - Baseline input-output data used as a reference.
#
# Process:
#   1. Computes ZX as the ratio of StockMat to E_CZ. This represents effective stock.
#   2. Creates a temporary matrix (temp) by replicating the first row of IOX_0 NN times,
#      scaled by 1.25.
#   3. Uses the global index matrix 'mat.key' to replace specific entries in ZX with
#      corresponding values from temp.
#   4. Computes VX as the ratio of IOV_t to E_VA, representing effective inflow.
#   5. Sums the rows of OrderMat to form a single-row matrix OX.
#   6. Combines ZX, VX, and OX by stacking them and returns the minimum value for each
#      sector-region (column).
#
# Returns:
#   A vector (or 1-row matrix) with length equal to RR*NN representing production output.
###############################################################################
Production = function(StockMat, E_CZ, IOV_t, E_VA, OrderMat, IOX_0)
{
    # Compute effective stock by dividing StockMat by efficiency E_CZ.
    ZX = StockMat / E_CZ   # ZX: Matrix with dimensions (NN) x (RR*NN)
    
    # Create a temporary matrix: replicate the first row of IOX_0 NN times and scale by 1.25.
    temp = 1.25 * IOX_0[rep(1, NN), ]
    
    # Replace elements in ZX (identified by the global index 'mat.key') with the corresponding
    # values from 'temp'. This injects specific baseline values into ZX.
    ZX[mat.key] = temp[mat.key]
    
    # Compute effective inflows/outflows by dividing IOV_t by efficiency E_VA.
    VX = IOV_t / E_VA   # VX: Matrix with dimensions (UU) x (RR*NN)
    
    # Sum orders across each row of OrderMat and transpose the result,
    # yielding a 1 x (RR*NN) vector.
    OX = t(apply(OrderMat, 1, sum))
    
    # Combine ZX, VX, and OX into a single matrix by row-binding and take the minimum
    # across each column. This ensures production is limited by the tightest constraint.
    return(apply(rbind(ZX, VX, OX), 2, min))
}

###############################################################################
# Production_max: Calculates the maximum possible production output based on
# available stock and inflows without the order constraint.
#
# Parameters: Same as Production.
#
# Process:
#   1. Computes ZX similarly to Production.
#   2. Replaces specific entries in ZX using the temporary matrix 'temp' and index 'mat.key'.
#   3. Computes VX as the ratio of IOV_t to E_VA.
#   4. Returns the minimum of ZX and VX (ignoring order constraints from OrderMat).
#
# Returns:
#   A vector (or 1-row matrix) representing the maximum possible production output.
###############################################################################
Production_max = function(StockMat, E_CZ, IOV_t, E_VA, IOX_0)
{
    # Compute effective stock.
    ZX = StockMat / E_CZ   # ZX: Matrix with dimensions (NN) x (RR*NN)
    
    # Create temporary matrix from the first row of IOX_0 (replicated NN times) scaled by 1.25.
    temp = 1.25 * IOX_0[rep(1, NN), ]
    
    # Replace entries in ZX as indicated by mat.key.
    ZX[mat.key] = temp[mat.key]
    
    # Compute effective inflow.
    VX = IOV_t / E_VA   # VX: Matrix with dimensions (UU) x (RR*NN)
    
    # Return the minimum of ZX and VX for each sector-region.
    return(apply(rbind(ZX, VX), 2, min))
}

###############################################################################
# OverProdSignFun: Determines the sign and magnitude of overproduction
# relative to the baseline, comparing actual orders to available capacity.
#
# Parameters:
#   StockMat  - Matrix representing the available stock.
#   E_Z       - Efficiency factor for the stock (similar to E_CZ but for this function).
#   IOV_t     - Matrix of inflows/outflows at time t.
#   E_V       - Efficiency factor for the inflows/outflows (analogous to E_VA).
#   OrderMat  - Matrix representing orders.
#   IOX_0     - Baseline input-output data used to scale the differences.
#
# Process:
#   1. Initializes a Result matrix of zeros to store the overproduction sign and magnitude.
#   2. Computes ZX as StockMat/E_Z (adjusted stock).
#   3. Computes VX as IOV_t/E_V (adjusted inflows).
#   4. Computes OX as the sum of orders for each sector-region.
#   5. Combines ZX and OX into a new matrix ZOX.
#   6. Compares VX (replicated across rows) with ZOX:
#       - If VX is less than ZOX for all rows in a column, that column gets a +1 sign.
#       - If any entry of VX exceeds ZOX, that column gets a -1 sign.
#   7. Computes a scaling factor based on the absolute difference between OX and VX,
#      normalized by IOX_0.
#   8. Returns the product of the sign matrix and the scaling factor.
#
# Returns:
#   A matrix of dimensions (UU) x (RR*NN) indicating the signed magnitude of overproduction.
###############################################################################
OverProdSignFun = function(StockMat, E_Z, IOV_t, E_V, OrderMat, IOX_0)
{
    # Initialize the result matrix with zeros.
    Result = matrix(0, nrow = UU, ncol = RR*NN)
    
    # Compute the adjusted stock level.
    ZX = StockMat / E_Z   # ZX: (RR*NN) x (RR*NN)
    
    # Compute the adjusted inflows.
    VX = IOV_t / E_V      # VX: (UU) x (RR*NN)
    
    # Sum the orders for each sector-region and transpose to get a 1 x (RR*NN) vector.
    OX = t(apply(OrderMat, 1, sum))
    
    # Combine ZX and OX into one matrix.
    ZOX = rbind(ZX, OX)
    
    # Create a logical matrix where each element indicates if the corresponding 
    # element of VX (replicated to the same dimensions as ZOX) is less than ZOX.
    V1U = VX[rep(1, nrow(ZOX)), ] < ZOX
    # If all values in a column are TRUE, then set that column's sign to +1.
    Result[1, which(apply(V1U, 2, sum) == nrow(ZOX))] = 1
    
    # Similarly, create a logical matrix for when VX exceeds ZOX.
    V1D = VX[rep(1, nrow(ZOX)), ] > ZOX
    # If any value in a column is TRUE, then set that column's sign to -1.
    Result[1, which(apply(V1D, 2, sum) > 0)] = -1
    
    # Compute the absolute difference between orders and inflows,
    # then normalize by IOX_0 (scaling factor from the baseline).
    temp = abs(OX - VX) / IOX_0
    
    # Return the element-wise product of the Result sign matrix and the scaling factor.
    return(Result * temp)
}