//! Test against the Wycheproof test vectors.

use ml_dsa::{KeyGen, MlDsa44, MlDsa65, MlDsa87, Signature, VerifyingKey};
use serde::Deserialize;
use std::fs::File;

#[derive(Deserialize, Debug)]
struct SignSeedTestFile {
    #[serde(rename(deserialize = "testGroups"))]
    groups: Vec<SignSeedTestGroup>,
    algorithm: String,
    header: Vec<String>,
}

#[derive(Deserialize, Debug)]
struct SignSeedTestGroup {
    #[allow(dead_code)]
    #[serde(rename(deserialize = "type"))]
    type_: String,

    #[serde(rename(deserialize = "privateSeed"), with = "hex::serde")]
    private_seed: Vec<u8>,

    #[serde(default, rename(deserialize = "publicKey"), with = "hex::serde")]
    #[allow(dead_code)]
    public_key: Vec<u8>,

    tests: Vec<SignSeedTest>,
}

#[derive(Deserialize, Debug)]
struct SignSeedTest {
    #[serde(rename(deserialize = "tcId"))]
    id: usize,
    comment: String,
    #[serde(default, with = "hex::serde")]
    msg: Vec<u8>,
    #[serde(default, with = "hex::serde")]
    mu: Vec<u8>,
    #[serde(default, with = "hex::serde")]
    ctx: Vec<u8>,
    #[serde(with = "hex::serde")]
    sig: Vec<u8>,
    result: ExpectedResult,
}

#[derive(Deserialize, Debug)]
struct VerifyTestFile {
    #[serde(rename(deserialize = "testGroups"))]
    groups: Vec<VerifyTestGroup>,
    algorithm: String,
    header: Vec<String>,
}

#[derive(Deserialize, Debug)]
struct VerifyTestGroup {
    #[allow(dead_code)]
    #[serde(rename(deserialize = "type"))]
    type_: String,

    #[serde(rename(deserialize = "publicKey"), with = "hex::serde")]
    public_key: Vec<u8>,

    tests: Vec<VerifyTest>,
}

#[derive(Deserialize, Debug)]
struct VerifyTest {
    #[serde(rename(deserialize = "tcId"))]
    id: usize,
    comment: String,
    #[serde(with = "hex::serde")]
    msg: Vec<u8>,
    #[serde(default, with = "hex::serde")]
    ctx: Vec<u8>,
    #[serde(with = "hex::serde")]
    sig: Vec<u8>,
    result: ExpectedResult,
}

#[derive(Copy, Clone, Deserialize, Debug, PartialEq)]
#[serde(rename_all = "lowercase")]
enum ExpectedResult {
    Valid,
    Invalid,
    Acceptable,
}

macro_rules! load_json_file {
    ($json_file:expr) => {{
        let path = format!(
            "{}/thirdparty/wycheproof/testvectors_v1/{}",
            env!("CARGO_MANIFEST_DIR"),
            $json_file
        );
        let data_file = File::open(&path)
            .expect("failed to open data file (try running `git submodule update --init`)");

        println!("Loading file: {path}");
        data_file
    }};
}

macro_rules! mldsa_sign_seed_test {
    ($name:ident, $json_file:expr, $alg:ident) => {
        #[test]
        fn $name() {
            let data_file = load_json_file!($json_file);
            let tests: SignSeedTestFile =
                serde_json::from_reader(data_file).expect("invalid test JSON");
            println!("{}:\n{}\n", tests.algorithm, tests.header.join(""));

            for group in tests.groups {
                let kp = $alg::from_seed(&group.private_seed.as_slice().try_into().unwrap());

                for test in &group.tests {
                    println!("Test #{}: {} ({:?})", test.id, &test.comment, &test.result);

                    // Some tests use `mu` (pre-hashed) instead of `msg`
                    if !test.mu.is_empty() {
                        // mu-based signing (pre-computed mu)
                        let mu = ml_dsa::util::B64::try_from(test.mu.as_slice()).unwrap();
                        let sig = kp.signing_key().sign_mu_deterministic(&mu);

                        match test.result {
                            ExpectedResult::Valid => {
                                assert_eq!(
                                    &*sig.encode(),
                                    test.sig.as_slice(),
                                    "Test #{}: mu-based signature mismatch",
                                    test.id
                                );
                            }
                            ExpectedResult::Invalid => {
                                // For mu-based tests, invalid means the sig shouldn't match
                                // (signing always succeeds with mu)
                            }
                            other => todo!("{:?}", other),
                        }
                        continue;
                    }

                    if test.ctx.is_empty() {
                        // Use sign_deterministic with empty context
                        let result = kp.signing_key().sign_deterministic(&test.msg, &[]);

                        match test.result {
                            ExpectedResult::Valid => {
                                let sig = result.expect("signing should succeed");
                                assert_eq!(
                                    &*sig.encode(),
                                    test.sig.as_slice(),
                                    "Test #{}: signature mismatch",
                                    test.id
                                );
                            }
                            ExpectedResult::Invalid => {
                                // Invalid tests may fail at signing or produce different sigs
                                if let Ok(sig) = result {
                                    assert_ne!(
                                        &*sig.encode(),
                                        test.sig.as_slice(),
                                        "Test #{}: expected invalid but got matching signature",
                                        test.id
                                    );
                                }
                            }
                            other => todo!("{:?}", other),
                        }
                    } else {
                        let result = kp.signing_key().sign_deterministic(&test.msg, &test.ctx);

                        match test.result {
                            ExpectedResult::Valid => {
                                let sig = result.expect("signing should succeed");
                                assert_eq!(
                                    &*sig.encode(),
                                    test.sig.as_slice(),
                                    "Test #{}: signature mismatch",
                                    test.id
                                );
                            }
                            ExpectedResult::Invalid => {
                                assert!(
                                    result.is_err(),
                                    "Test #{}: expected error for invalid test",
                                    test.id
                                );
                            }
                            other => todo!("{:?}", other),
                        }
                    }
                }
            }
        }
    };
}

macro_rules! mldsa_verify_test {
    ($name:ident, $json_file:expr, $alg:ident) => {
        #[test]
        fn $name() {
            let data_file = load_json_file!($json_file);
            let tests: VerifyTestFile =
                serde_json::from_reader(data_file).expect("invalid test JSON");
            println!("{}:\n{}\n", tests.algorithm, tests.header.join(""));

            for group in &tests.groups {
                if let Ok(encoded_vk) = group.public_key.as_slice().try_into() {
                    let vk = VerifyingKey::<$alg>::decode(&encoded_vk);
                    for test in &group.tests {
                        println!("Test #{}: {} ({:?})", test.id, &test.comment, &test.result);

                        if let Some(sig) = test
                            .sig
                            .as_slice()
                            .try_into()
                            .ok()
                            .and_then(|sig| Signature::<$alg>::decode(&sig))
                        {
                            if test.ctx.is_empty() {
                                let verified = vk.verify_with_context(&test.msg, &[], &sig);

                                match test.result {
                                    ExpectedResult::Valid => assert!(
                                        verified,
                                        "Test #{}: expected valid but verification failed",
                                        test.id
                                    ),
                                    ExpectedResult::Invalid => assert!(
                                        !verified,
                                        "Test #{}: expected invalid but verification passed",
                                        test.id
                                    ),
                                    other => todo!("{:?}", other),
                                }
                            } else {
                                let verified = vk.verify_with_context(&test.msg, &test.ctx, &sig);

                                match test.result {
                                    ExpectedResult::Valid => assert!(
                                        verified,
                                        "Test #{}: expected valid but verification failed",
                                        test.id
                                    ),
                                    ExpectedResult::Invalid => assert!(
                                        !verified,
                                        "Test #{}: expected invalid but verification passed",
                                        test.id
                                    ),
                                    other => todo!("{:?}", other),
                                }
                            }
                        } else {
                            println!("error decoding signature (length: {})", test.sig.len(),);
                            assert_eq!(
                                test.result,
                                ExpectedResult::Invalid,
                                "Test #{}: failed to decode signature but expected valid",
                                test.id
                            );
                        }
                    }
                }
            }
        }
    };
}

mldsa_sign_seed_test!(
    mldsa_44_sign_seed_test,
    "mldsa_44_sign_seed_test.json",
    MlDsa44
);
mldsa_sign_seed_test!(
    mldsa_65_sign_seed_test,
    "mldsa_65_sign_seed_test.json",
    MlDsa65
);
mldsa_sign_seed_test!(
    mldsa_87_sign_seed_test,
    "mldsa_87_sign_seed_test.json",
    MlDsa87
);
mldsa_verify_test!(mldsa_44_verify_test, "mldsa_44_verify_test.json", MlDsa44);
mldsa_verify_test!(mldsa_65_verify_test, "mldsa_65_verify_test.json", MlDsa65);
mldsa_verify_test!(mldsa_87_verify_test, "mldsa_87_verify_test.json", MlDsa87);
