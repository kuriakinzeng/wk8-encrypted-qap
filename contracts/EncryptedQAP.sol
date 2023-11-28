// SPDX-License-Identifier: SEE LICENSE IN LICENSE
pragma solidity ^0.8.10;

contract EncryptedQAP {
    struct ECPoint {
        uint256 x;
        uint256 y;
    } 

    struct ECPoint2 {
        uint256[2] x;
        uint256[2] y;
    } 

    uint256 constant FIELD_MOD = 21888242871839275222246405745257275088696311157297823662689037894645226208583;

    function neg(ECPoint memory pt) internal pure returns (ECPoint memory) {
        if (pt.x == 0 && pt.y == 0)
            return pt;
        else
            return ECPoint(pt.x, (FIELD_MOD - pt.y) % FIELD_MOD);
    }

    function verify(ECPoint memory A1, ECPoint2 memory B2, ECPoint memory C1) public view returns (bool) {
        ECPoint2 memory G2 = ECPoint2(
            [   
                10857046999023057135944570762232829481370756359578518086990519993285655852781, 
                11559732032986387107991004021392285783925812861821192530917403151452391805634
            ],
            [
                8495653923123431417604973247489272438418190587263600148770280649306958101930, 
                4082367875863433681332203403145435568316851327593401208105741076214120093531
            ]
        );

        ECPoint memory negA1 = neg(A1);

        uint256[12] memory points = [
            negA1.x,
            negA1.y,
            B2.x[1],
            B2.x[0],
            B2.y[1],
            B2.y[0],
            C1.x,
            C1.y,
            G2.x[1],
            G2.x[0],
            G2.y[1],
            G2.y[0]
        ];

        bool success;

        assembly {
            success := staticcall(gas(), 0x08, points, mul(12, 0x20), points, 0x20)
        }
        
        return success;
    }
}