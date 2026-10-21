// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "transit/image.h"

namespace transit {

struct HeaderCard {
    enum class Type { String, Int, Float, Bool, Undefined };
    std::string key, raw_value, comment;
    Type type = Type::Undefined;
};

/// Backend-neutral reader. The pipeline only talks to this interface, so the
/// CFITSIO and libxisfconv backends can be swapped (and compared) per run.
class ImageReader {
public:
    virtual ~ImageReader() = default;
    virtual void open(const std::string& path) = 0;
    virtual int hdu_count() const = 0;
    virtual std::vector<HeaderCard> header(int hdu = 0) const = 0;
    virtual std::vector<long> dims(int hdu = 0) const = 0;
    /// Physical values: BSCALE/BZERO applied, converted to float.
    virtual Image pixels(int hdu = 0) const = 0;
    virtual std::string backend() const = 0;
};

/// "cfitsio" or "xisfconv"; throws if the backend was not compiled in.
std::unique_ptr<ImageReader> make_reader(std::string_view backend);
std::vector<std::string> available_backends();
/// True if separate readers of this backend may be used from different threads at the same time.
bool backend_thread_safe(std::string_view backend);

/// Write a 32-bit float primary image (CFITSIO builds only; throws otherwise).
void write_fits_float(const std::string& path, const Image& img, const std::vector<HeaderCard>& extra = {});

/// Convenience lookups on a parsed header.
std::optional<std::string> card_string(const std::vector<HeaderCard>& h, std::string_view key);
std::optional<double> card_double(const std::vector<HeaderCard>& h, std::string_view key);

}  // namespace transit
