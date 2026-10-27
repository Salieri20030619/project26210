#!/usr/bin/env bash
set -euo pipefail

if [[ "$#" -ne 4 ]]; then
  echo "usage: download_zenodo_segmented_v1.sh OUTPUT_DIR FILE_NAME EXPECTED_SIZE EXPECTED_MD5" >&2
  exit 2
fi

output_dir="$1"
file_name="$2"
expected_size="$3"
expected_md5="$4"
segments=8
url="https://zenodo.org/api/records/11286535/files/${file_name}/content"
part_dir="${output_dir}/.${file_name}.segments_v1"
assembled="${output_dir}/.${file_name}.assembled_v1"
final="${output_dir}/${file_name}"

mkdir -p "$output_dir" "$part_dir"
if [[ -f "$final" ]]; then
  observed_size="$(stat -c %s "$final")"
  observed_md5="$(md5sum "$final" | awk '{print $1}')"
  if [[ "$observed_size" == "$expected_size" && "$observed_md5" == "$expected_md5" ]]; then
    echo "already valid: $file_name"
    exit 0
  fi
  mv "$final" "${output_dir}/.${file_name}.serial_partial_v1"
fi

chunk=$(( (expected_size + segments - 1) / segments ))
download_segment() {
  local i="$1"
  start=$(( i * chunk ))
  end=$(( start + chunk - 1 ))
  if (( end >= expected_size )); then end=$(( expected_size - 1 )); fi
  if (( start > end )); then return 0; fi
  part="${part_dir}/part_${i}.bin"
  expected_part=$(( end - start + 1 ))
  if [[ -f "$part" ]] && (( $(stat -c %s "$part") > expected_part )); then
    echo "oversized segment: $part" >&2
    return 1
  fi
  while [[ ! -f "$part" ]] || (( $(stat -c %s "$part") < expected_part )); do
    observed=0
    if [[ -f "$part" ]]; then observed="$(stat -c %s "$part")"; fi
    resume_start=$(( start + observed ))
    tail="${part}.next"
    rm -f "$tail"
    set +e
    curl -sS -L --fail --retry 5 --retry-delay 2 --range "${resume_start}-${end}" "$url" -o "$tail"
    rc=$?
    set -e
    got=0
    if [[ -f "$tail" ]]; then got="$(stat -c %s "$tail")"; fi
    if (( got > expected_part - observed )); then
      echo "range response exceeded expected bytes: $part" >&2
      return 1
    fi
    if (( got > 0 )); then
      dd if="$tail" of="$part" oflag=append conv=notrunc status=none
    fi
    rm -f "$tail"
    if (( rc != 0 && got == 0 )); then sleep 2; fi
  done
  echo "segment PASS: $i bytes=$(stat -c %s "$part")"
}

pids=()
for ((i=0; i<segments; i++)); do
  download_segment "$i" &
  pids+=("$!")
done
for pid in "${pids[@]}"; do wait "$pid"; done

rm -f "$assembled"
for ((i=0; i<segments; i++)); do
  part="${part_dir}/part_${i}.bin"
  [[ -f "$part" ]] || continue
  dd if="$part" of="$assembled" oflag=append conv=notrunc status=none
done

observed_size="$(stat -c %s "$assembled")"
observed_md5="$(md5sum "$assembled" | awk '{print $1}')"
if [[ "$observed_size" != "$expected_size" || "$observed_md5" != "$expected_md5" ]]; then
  echo "validation failed: $file_name size=$observed_size md5=$observed_md5" >&2
  exit 1
fi
mv "$assembled" "$final"
echo "PASS $file_name size=$observed_size md5=$observed_md5"
